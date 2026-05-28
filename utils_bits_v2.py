# utils_bits.py
from __future__ import annotations

import hashlib
import os
import re
import zlib
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple


# ----------------------------
# Hash / IO helpers
# ----------------------------

def sha256_bytes(data: bytes) -> str:
    h = hashlib.sha256()
    h.update(data)
    return h.hexdigest()


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_dir(path: str) -> None:
    if path:
        os.makedirs(path, exist_ok=True)


def write_bytes(path: str, data: bytes) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "wb") as f:
        f.write(data)


def write_text(path: str, text: str) -> None:
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def read_bytes(path: str) -> bytes:
    with open(path, "rb") as f:
        return f.read()


# ----------------------------
# Bytes <-> bitstring
# ----------------------------

def bytes_to_bitstring(data: bytes) -> str:
    """Convert bytes to a '0'/'1' bitstring (MSB-first per byte)."""
    return "".join(f"{b:08b}" for b in data)


def bitstring_to_bytes(bitstr: str, pad_to_byte: bool = True) -> Tuple[bytes, int]:
    """
    Convert '0'/'1' bitstring (MSB-first) to bytes.
    Returns (bytes, pad_bits_added).
    """
    if bitstr is None:
        return b"", 0
    s = bitstr.strip()
    if s == "":
        return b"", 0
    if any(c not in "01" for c in s):
        raise ValueError("bitstr must contain only '0'/'1'")

    pad_bits = 0
    if (len(s) % 8) != 0:
        if not pad_to_byte:
            raise ValueError("bitstr length is not multiple of 8")
        pad_bits = 8 - (len(s) % 8)
        s = s + ("0" * pad_bits)

    out = bytearray()
    for i in range(0, len(s), 8):
        out.append(int(s[i:i+8], 2))
    return bytes(out), pad_bits


# ----------------------------
# Magic detection (headerless routing)
# ----------------------------

@dataclass
class MagicInfo:
    kind: str
    ext: str
    mime: str
    confidence: float
    note: str = ""


def _zip_kind_from_members(names: list[str]) -> Optional[MagicInfo]:
    """
    Heuristic: detect OOXML family (docx/pptx/xlsx) from ZIP members.
    """
    s = set(names)
    if "[Content_Types].xml" not in s:
        return None

    # OOXML folders
    if any(n.startswith("word/") for n in s):
        return MagicInfo(kind="docx", ext=".docx", mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document", confidence=0.95)
    if any(n.startswith("ppt/") for n in s):
        return MagicInfo(kind="pptx", ext=".pptx", mime="application/vnd.openxmlformats-officedocument.presentationml.presentation", confidence=0.95)
    if any(n.startswith("xl/") for n in s):
        return MagicInfo(kind="xlsx", ext=".xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", confidence=0.95)

    # EPUB is also a ZIP container
    if "mimetype" in s and any(n.startswith("META-INF/") for n in s):
        return MagicInfo(kind="epub", ext=".epub", mime="application/epub+zip", confidence=0.85)

    return None


def detect_magic(data: bytes) -> Optional[MagicInfo]:
    """
    Detect common container/codec signatures so decoder can route without a custom header.
    This is intentionally conservative; extend as needed.
    """
    if data is None or len(data) < 4:
        return None

    head = data[:64]

    # PDF
    if head.startswith(b"%PDF-"):
        return MagicInfo(kind="pdf", ext=".pdf", mime="application/pdf", confidence=0.97)

    # GZIP
    if head.startswith(b"\x1F\x8B"):
        return MagicInfo(kind="gzip", ext=".gz", mime="application/gzip", confidence=0.99)

    # PNG
    if head.startswith(b"\x89PNG\r\n\x1a\n"):
        return MagicInfo(kind="png", ext=".png", mime="image/png", confidence=0.99)

    # JPEG
    if head.startswith(b"\xFF\xD8\xFF"):
        return MagicInfo(kind="jpeg", ext=".jpg", mime="image/jpeg", confidence=0.98)

    # GIF
    if head.startswith(b"GIF87a") or head.startswith(b"GIF89a"):
        return MagicInfo(kind="gif", ext=".gif", mime="image/gif", confidence=0.98)

    # BMP
    if head.startswith(b"BM"):
        return MagicInfo(kind="bmp", ext=".bmp", mime="image/bmp", confidence=0.85)

    # TIFF
    if head.startswith(b"II*\x00") or head.startswith(b"MM\x00*"):
        return MagicInfo(kind="tiff", ext=".tif", mime="image/tiff", confidence=0.90)

    # WebP: RIFF....WEBP
    if head.startswith(b"RIFF") and len(head) >= 12 and head[8:12] == b"WEBP":
        return MagicInfo(kind="webp", ext=".webp", mime="image/webp", confidence=0.98)

    # WAV / AVI: RIFF....WAVE / RIFF....AVI
    if head.startswith(b"RIFF") and len(head) >= 12:
        tag = head[8:12]
        if tag == b"WAVE":
            return MagicInfo(kind="wav", ext=".wav", mime="audio/wav", confidence=0.98)
        if tag == b"AVI ":
            return MagicInfo(kind="avi", ext=".avi", mime="video/x-msvideo", confidence=0.85)

    # MP3
    if head.startswith(b"ID3"):
        return MagicInfo(kind="mp3", ext=".mp3", mime="audio/mpeg", confidence=0.85, note="ID3")
    if len(head) >= 2 and head[0] == 0xFF and (head[1] & 0xE0) == 0xE0:
        return MagicInfo(kind="mp3", ext=".mp3", mime="audio/mpeg", confidence=0.65, note="frame_sync")

    # Ogg container
    if head.startswith(b"OggS"):
        if b"OpusHead" in data[:4096]:
            return MagicInfo(kind="opus_ogg", ext=".ogg", mime="audio/ogg", confidence=0.95, note="OpusHead")
        return MagicInfo(kind="ogg", ext=".ogg", mime="application/ogg", confidence=0.90)

    # FLAC
    if head.startswith(b"fLaC"):
        return MagicInfo(kind="flac", ext=".flac", mime="audio/flac", confidence=0.99)

    # MP4/ISOBMFF: 'ftyp' at offset 4
    if len(head) >= 12 and head[4:8] == b"ftyp":
        return MagicInfo(kind="mp4", ext=".mp4", mime="video/mp4", confidence=0.92)

    # Matroska / WebM (EBML)
    if head.startswith(b"\x1A\x45\xDF\xA3"):
        return MagicInfo(kind="mkv_webm", ext=".mkv", mime="video/x-matroska", confidence=0.60, note="EBML")

    # ZIP: PK...
    if head.startswith(b"PK\x03\x04") or head.startswith(b"PK\x05\x06") or head.startswith(b"PK\x07\x08"):
        # Try to distinguish docx/pptx/xlsx vs generic zip using members (no custom header needed).
        try:
            import io, zipfile
            with zipfile.ZipFile(io.BytesIO(data), "r") as zf:
                names = zf.namelist()
            mk = _zip_kind_from_members(names)
            if mk is not None:
                return mk
        except Exception:
            pass
        return MagicInfo(kind="zip", ext=".zip", mime="application/zip", confidence=0.99)

    # Plain text heuristic (UTF-8)
    try:
        s = data[:4096].decode("utf-8")
        printable = sum(ch.isprintable() or ch in "\r\n\t" for ch in s)
        ratio = printable / max(1, len(s))
        if ratio > 0.95:
            return MagicInfo(kind="text", ext=".txt", mime="text/plain", confidence=0.60, note="heuristic")
    except Exception:
        pass

    return None


# ----------------------------
# ZLIB wrapper (headerless framing)
# ----------------------------

def _looks_already_compressed(m: Optional[MagicInfo], data: bytes) -> bool:
    if m and m.kind in {"zip","docx","pptx","xlsx","epub","png","jpeg","webp","gzip","flac","opus_ogg","ogg","mp4","pdf","mp3","wav","avi","mkv_webm","gif","bmp","tiff"}:
        return True
    if len(data) >= 2048:
        sample = data[:2048]
        uniq = len(set(sample))
        if uniq > 200:
            return True
    return False


def zlib_wrap(inner_bytes: bytes, policy: str = "auto") -> Tuple[bytes, Dict[str, Any]]:
    """
    Return zlib stream bytes S that is self-terminating + self-checking.
    policy:
      - "stored": force level=0 (prefer STORED blocks)
      - "compress": force compression (level=6)
      - "auto": stored if already compressed, else compress
    """
    m = detect_magic(inner_bytes)
    if policy == "stored":
        level = 0
        why = "forced_stored"
    elif policy == "compress":
        level = 6
        why = "forced_compress"
    else:
        if _looks_already_compressed(m, inner_bytes):
            level = 0
            why = "auto_stored_already_compressed"
        else:
            level = 6
            why = "auto_compress"

    z = zlib.compress(inner_bytes, level=level)  # zlib header + adler32
    meta = {
        "zlib_level": level,
        "policy": policy,
        "decision": why,
        "inner_magic": (m.kind if m else None),
    }
    return z, meta


def zlib_inflate_until_eof(buffer: bytes) -> Tuple[bytes, Dict[str, Any]]:
    """
    Inflate ONLY the first valid zlib stream in buffer.
    Trailing bytes (from DNA decode padding, or other residue) are ignored via unused_data.
    Integrity is validated by zlib checksum; if corrupt => error.
    """
    dco = zlib.decompressobj(wbits=zlib.MAX_WBITS)
    out = bytearray()
    try:
        out.extend(dco.decompress(buffer))
        while dco.unconsumed_tail and not dco.eof:
            out.extend(dco.decompress(dco.unconsumed_tail))
        out.extend(dco.flush())
        meta = {
            "eof": bool(dco.eof),
            "unused_tail_len": len(dco.unused_data),
            "unconsumed_tail_len": len(dco.unconsumed_tail),
        }
        return bytes(out), meta
    except zlib.error as e:
        return b"", {"eof": False, "error": str(e), "unused_tail_len": 0, "unconsumed_tail_len": 0}


# ----------------------------
# Filename sanitization
# ----------------------------

_INVALID_CHARS = re.compile(r"[^A-Za-z0-9._-]+")


def safe_basename(name: str, fallback: str = "file.bin") -> str:
    base = os.path.basename(name) if name else fallback
    base = base.strip().replace(" ", "_")
    base = _INVALID_CHARS.sub("_", base)
    if not base:
        base = fallback
    return base
