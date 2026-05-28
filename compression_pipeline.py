# compression_pipeline.py
# Headerless self-describing compression benchmark
# ------------------------------------------------
# New policy:
#   - No Lossy/Lossless UI dependency.
#   - Compression mode benchmarks all available valid compressors.
#   - Only self-describing outputs are accepted, based on detect_magic().
#   - The smallest byte stream is selected for DNA encoding.
#
# Requirements from existing project:
#   - utils_bits_v2.py: detect_magic, safe_basename
#   - ui_helpers.py: fmt_bytes, get_domain
#   - config.py: SELF_DESCRIBING_KINDS
#
# Optional:
#   - Pillow for image codecs
#   - ffmpeg on PATH for audio/video codecs

from __future__ import annotations

import bz2
import gzip
import hashlib
import io
import lzma
import os
import subprocess
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

try:
    from PIL import Image
except Exception:
    Image = None

from utils_bits_v2 import detect_magic, safe_basename
from config import SELF_DESCRIBING_KINDS
from ui_helpers import fmt_bytes, get_domain


# ============================================================
# Data model
# ============================================================

@dataclass
class CompressionCandidate:
    rank: int
    method: str
    data: bytes
    ext: str
    kind: str
    mime: str
    lossy: bool
    size_bytes: int
    compression_ratio: float
    saving_pct: float
    estimated_dna_nt: int
    note: str = ""

    def public_row(self) -> dict:
        return {
            "Rank": self.rank,
            "Method": self.method,
            
            "File extension": self.ext,
            
            "Size": fmt_bytes(self.size_bytes),
            "Ratio": f"{self.compression_ratio:.2f}",
           
            "Estimate DNA (nt)": f"{self.estimated_dna_nt:,}"
        }


def _candidate_from_bytes(
    method: str,
    data: bytes,
    original_size: int,
    lossy: bool,
    note: str = "",
) -> Optional[CompressionCandidate]:
    """Create a candidate only if bytes are self-describing by magic signature."""
    if not data:
        return None

    m = detect_magic(data)
    if not m:
        return None

    if m.kind not in SELF_DESCRIBING_KINDS:
        return None

    size = len(data)
    ratio = (original_size / size) if size else 0.0
    saving = (1.0 - (size / original_size)) * 100.0 if original_size else 0.0

    return CompressionCandidate(
        rank=0,
        method=method,
        data=data,
        ext=m.ext,
        kind=m.kind,
        mime=m.mime,
        lossy=lossy,
        size_bytes=size,
        compression_ratio=ratio,
        saving_pct=saving,
        estimated_dna_nt=size * 4,  # exact for simple 2-bit mapping; estimate for rule-based mapping
        note=note or getattr(m, "note", ""),
    )


# ============================================================
# Generic self-describing containers
# ============================================================

def zip_store_bytes(name: str, data: bytes) -> bytes:
    """ZIP container with no compression; useful to force self-describing bytes."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_STORED) as zf:
        zf.writestr(safe_basename(name), data)
    return buf.getvalue()


def zip_deflate_bytes(name: str, data: bytes, level: int = 6) -> bytes:
    """ZIP container with DEFLATE compression."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=int(level)) as zf:
        zf.writestr(safe_basename(name), data)
    return buf.getvalue()


def generic_container_candidates(raw: bytes, original_name: str, original_size: int) -> List[CompressionCandidate]:
    """
    Generic lossless self-describing candidates.

    Used for:
      - text/json/csv
      - binary/unknown
      - fallback for image/audio/video/document/archive
    """
    out: List[CompressionCandidate] = []

    # gzip
    for lvl in (1, 6, 9):
        try:
            b = gzip.compress(raw, compresslevel=int(lvl))
            c = _candidate_from_bytes(f"gzip_lvl{lvl}", b, original_size, lossy=False, note="generic lossless")
            if c:
                out.append(c)
        except Exception:
            pass

    # bz2
    for lvl in (1, 6, 9):
        try:
            b = bz2.compress(raw, compresslevel=int(lvl))
            c = _candidate_from_bytes(f"bz2_lvl{lvl}", b, original_size, lossy=False, note="generic lossless")
            if c:
                out.append(c)
        except Exception:
            pass

    # xz
    for preset in (0, 6, 9):
        try:
            b = lzma.compress(raw, format=lzma.FORMAT_XZ, preset=int(preset))
            c = _candidate_from_bytes(f"xz_p{preset}", b, original_size, lossy=False, note="generic lossless")
            if c:
                out.append(c)
        except Exception:
            pass

    # zip store
    try:
        b = zip_store_bytes(original_name, raw)
        c = _candidate_from_bytes("zip_store", b, original_size, lossy=False, note="self-describing wrapper")
        if c:
            out.append(c)
    except Exception:
        pass

    # zip deflate
    for lvl in (1, 6, 9):
        try:
            b = zip_deflate_bytes(original_name, raw, level=int(lvl))
            c = _candidate_from_bytes(f"zip_deflate_lvl{lvl}", b, original_size, lossy=False, note="generic lossless")
            if c:
                out.append(c)
        except Exception:
            pass

    return out


# ============================================================
# Image candidates
# ============================================================

def image_candidates(raw: bytes, original_name: str, original_size: int, quality_mode: str = "Compression") -> List[CompressionCandidate]:
    """
    Image-native candidates.

    In the new pipeline, Compression mode always tries both:
      - lossless image candidates: PNG, WebP lossless
      - lossy image candidates: WebP q50/q60/q70/q80/q90, JPEG q50/q60/q70/q80/q90

    quality_mode is kept only for backward compatibility and is not used.
    """
    out: List[CompressionCandidate] = []
    if Image is None:
        return out

    try:
        img = Image.open(io.BytesIO(raw))
        img.load()
    except Exception:
        return out

    # PNG lossless
    try:
        im = img
        if im.mode not in ("RGB", "RGBA", "L", "P"):
            im = im.convert("RGBA")
        buf = io.BytesIO()
        im.save(buf, format="PNG", optimize=True)
        c = _candidate_from_bytes("png_lossless", buf.getvalue(), original_size, lossy=False, note="image-native lossless")
        if c:
            out.append(c)
    except Exception:
        pass

    # WebP lossless
    try:
        buf = io.BytesIO()
        img.save(buf, format="WEBP", lossless=True, quality=100)
        c = _candidate_from_bytes("webp_lossless", buf.getvalue(), original_size, lossy=False, note="image-native lossless")
        if c:
            out.append(c)
    except Exception:
        pass

    # WebP / JPEG lossy candidates
    for q in (50, 60, 70, 80, 90):
        # WebP lossy
        try:
            buf = io.BytesIO()
            img.save(buf, format="WEBP", quality=int(q), lossless=False)
            c = _candidate_from_bytes(f"webp_q{q}", buf.getvalue(), original_size, lossy=True, note="image-native lossy")
            if c:
                out.append(c)
        except Exception:
            pass

        # JPEG lossy
        try:
            im = img
            if im.mode != "RGB":
                im = im.convert("RGB")
            buf = io.BytesIO()
            im.save(buf, format="JPEG", quality=int(q), optimize=True)
            c = _candidate_from_bytes(f"jpeg_q{q}", buf.getvalue(), original_size, lossy=True, note="image-native lossy")
            if c:
                out.append(c)
        except Exception:
            pass

    return out


# ============================================================
# ffmpeg helpers
# ============================================================

def has_ffmpeg() -> bool:
    """Return True if ffmpeg is available on PATH."""
    try:
        p = subprocess.run(
            ["ffmpeg", "-version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=5,
        )
        return p.returncode == 0
    except Exception:
        return False


def run_ffmpeg_transcode(input_path: str, output_ext: str, args: List[str], timeout_sec: int = 120) -> Optional[bytes]:
    """
    Run ffmpeg and return output bytes.
    Returns None if ffmpeg is unavailable, codec fails, or timeout occurs.
    """
    try:
        with tempfile.TemporaryDirectory() as td:
            out_path = Path(td) / f"out{output_ext}"
            cmd = [
                "ffmpeg",
                "-y",
                "-hide_banner",
                "-loglevel",
                "error",
                "-i",
                input_path,
            ] + args + [str(out_path)]

            p = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=int(timeout_sec),
            )

            if p.returncode != 0 or not out_path.exists():
                return None

            return out_path.read_bytes()
    except Exception:
        return None


# ============================================================
# Audio candidates
# ============================================================

def audio_compression_candidates(
    input_path: str,
    raw: bytes,
    original_size: int,
    quality_mode: str = "Compression",
) -> List[CompressionCandidate]:
    """
    Audio-native candidates.

    New pipeline: always try audio-native candidates in Compression mode.

    Requires ffmpeg:
      - FLAC lossless
      - OGG/Opus lossy
      - MP3 lossy
      - AAC/M4A lossy
    """
    out: List[CompressionCandidate] = []

    if not has_ffmpeg():
        return out

    # FLAC lossless
    flac = run_ffmpeg_transcode(input_path, ".flac", ["-vn", "-c:a", "flac"])
    if flac:
        c = _candidate_from_bytes("flac_lossless", flac, original_size, lossy=False, note="audio-native lossless")
        if c:
            out.append(c)

    # OGG/Opus lossy
    for br in (32, 64, 96, 128):
        ogg = run_ffmpeg_transcode(
            input_path,
            ".ogg",
            ["-vn", "-c:a", "libopus", "-b:a", f"{br}k"],
        )
        if ogg:
            c = _candidate_from_bytes(f"opus_ogg_{br}k", ogg, original_size, lossy=True, note="audio-native lossy")
            if c:
                out.append(c)

    # MP3 lossy
    for br in (96, 128, 192):
        mp3 = run_ffmpeg_transcode(
            input_path,
            ".mp3",
            ["-vn", "-c:a", "libmp3lame", "-b:a", f"{br}k"],
        )
        if mp3:
            c = _candidate_from_bytes(f"mp3_{br}k", mp3, original_size, lossy=True, note="audio-native lossy")
            if c:
                out.append(c)

    # AAC/M4A lossy. detect_magic usually identifies this as mp4/ftyp.
    for br in (96, 128):
        aac = run_ffmpeg_transcode(
            input_path,
            ".m4a",
            ["-vn", "-c:a", "aac", "-b:a", f"{br}k"],
        )
        if aac:
            c = _candidate_from_bytes(f"aac_m4a_{br}k", aac, original_size, lossy=True, note="audio-native lossy")
            if c:
                out.append(c)

    return out


# ============================================================
# Video candidates
# ============================================================

def video_compression_candidates(
    input_path: str,
    raw: bytes,
    original_size: int,
    quality_mode: str = "Compression",
) -> List[CompressionCandidate]:
    """
    Video-native candidates.

    New pipeline: always try video-native candidates in Compression mode.

    Requires ffmpeg:
      - H.264 MP4 CRF 28/32
      - VP9 WebM CRF 32
    """
    out: List[CompressionCandidate] = []

    if not has_ffmpeg():
        return out

    # H.264 MP4
    for crf in (28, 32):
        mp4 = run_ffmpeg_transcode(
            input_path,
            ".mp4",
            [
                "-c:v", "libx264",
                "-preset", "veryfast",
                "-crf", str(crf),
                "-pix_fmt", "yuv420p",
                "-movflags", "+faststart",
                "-c:a", "aac",
                "-b:a", "128k",
            ],
            timeout_sec=180,
        )
        if mp4:
            c = _candidate_from_bytes(f"h264_mp4_crf{crf}", mp4, original_size, lossy=True, note="video-native lossy")
            if c:
                out.append(c)

    # VP9 WebM
    webm = run_ffmpeg_transcode(
        input_path,
        ".webm",
        [
            "-c:v", "libvpx-vp9",
            "-crf", "32",
            "-b:v", "0",
            "-row-mt", "1",
            "-cpu-used", "4",
            "-pix_fmt", "yuv420p",
            "-c:a", "libopus",
            "-b:a", "96k",
        ],
        timeout_sec=180,
    )
    if webm:
        c = _candidate_from_bytes("vp9_webm_crf32", webm, original_size, lossy=True, note="video-native lossy")
        if c:
            out.append(c)

    return out


# ============================================================
# Benchmark orchestration
# ============================================================

def _deduplicate_candidates(candidates: List[CompressionCandidate]) -> List[CompressionCandidate]:
    """Deduplicate by method, size, and sha1 prefix."""
    seen = set()
    unique: List[CompressionCandidate] = []

    for c in candidates:
        key = (c.method, c.size_bytes, hashlib.sha1(c.data[:4096]).hexdigest())
        if key in seen:
            continue
        seen.add(key)
        unique.append(c)

    return unique


def run_compression_benchmark(
    input_path: str,
    raw: bytes,
    quality_mode: str = "Compression",
) -> Tuple[CompressionCandidate, List[CompressionCandidate]]:
    """
    Build and rank all valid self-describing compression candidates.

    The smallest candidate by byte size is selected.

    quality_mode is retained for backward compatibility but should no longer
    control compressor inclusion. Compression mode tries all available candidates.
    """
    original_name = os.path.basename(input_path) or "input.bin"
    original_size = len(raw)
    domain = get_domain(input_path, raw)

    candidates: List[CompressionCandidate] = []

    # Keep original if it is already self-describing.
    keep = _candidate_from_bytes("keep_original", raw, original_size, lossy=False, note="native/self-describing input")
    if keep:
        candidates.append(keep)

    # Domain-specific candidates.
    if domain == "image":
        candidates.extend(image_candidates(raw, original_name, original_size, quality_mode=quality_mode))
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    elif domain == "audio":
        candidates.extend(audio_compression_candidates(input_path, raw, original_size, quality_mode=quality_mode))
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    elif domain == "video":
        candidates.extend(video_compression_candidates(input_path, raw, original_size, quality_mode=quality_mode))
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    elif domain in {"text", "other", "unknown", "binary"}:
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    elif domain in {"archive", "document"}:
        # Already-compressed/self-describing files often select keep_original,
        # but generic containers are included for comparison.
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    else:
        candidates.extend(generic_container_candidates(raw, original_name, original_size))

    unique = _deduplicate_candidates(candidates)

    # Last-resort ZIP store: ensures one self-describing payload even for raw/unknown bytes.
    if not unique:
        b = zip_store_bytes(original_name, raw)
        c = _candidate_from_bytes("zip_store_fallback", b, original_size, lossy=False, note="fallback wrapper")
        if not c:
            raise RuntimeError("Could not create any self-describing compression candidate.")
        unique = [c]

    # Rank by smallest stored output; this is the byte stream that becomes DNA.
    unique.sort(key=lambda x: (x.size_bytes, x.method))

    for i, c in enumerate(unique, start=1):
        c.rank = i

    return unique[0], unique
