from __future__ import annotations

import io
import os
import hashlib
import uuid
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import streamlit as st

try:
    from PIL import Image
except Exception:
    Image = None

import dna_codec
from utils_bits_v2 import detect_magic, safe_basename

try:
    from compressors_v2 import detect_domain
except Exception:
    detect_domain = None

from config import (
    WORK_ROOT,
    IMAGE_KINDS,
    AUDIO_KINDS,
    VIDEO_KINDS,
    IMAGE_PREVIEW_USE_CONTAINER_WIDTH,
    IMAGE_PREVIEW_WIDTH,
    TEXT_PREVIEW_HEIGHT,
)

_PREVIEW_FILE_CALL_COUNTER = 0


def fmt_bytes(n: Optional[int]) -> str:
    if n is None:
        return "—"
    try:
        x = float(n)
    except Exception:
        return "—"
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if x < 1024.0 or unit == "TB":
            return f"{int(x)} B" if unit == "B" else f"{x:.2f} {unit}"
        x /= 1024.0
    return f"{x:.2f} TB"

def step_header(number: int, title: str) -> None:
    st.markdown(
        f"""
<div class="step-heading">
  <span class="step-badge">{number}</span>
  <span class="step-title">{title}</span>
</div>
""",
        unsafe_allow_html=True,
    )

def safe_write_bytes(path: str | Path, data: bytes) -> str:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data)
    return str(p)

def save_upload(uploaded_file) -> Tuple[str, bytes]:
    run_dir = WORK_ROOT / "uploads" / uuid.uuid4().hex
    run_dir.mkdir(parents=True, exist_ok=True)
    name = safe_basename(uploaded_file.name or "upload.bin")
    data = uploaded_file.getvalue()
    path = run_dir / name
    path.write_bytes(data)
    return str(path), data

def magic_dict(data: bytes) -> Dict[str, Any]:
    m = detect_magic(data)
    if not m:
        return {"kind": "unknown", "ext": ".bin", "mime": "application/octet-stream", "confidence": 0.0, "note": ""}
    return {
        "kind": m.kind,
        "ext": m.ext,
        "mime": m.mime,
        "confidence": m.confidence,
        "note": getattr(m, "note", ""),
    }

def get_domain(path: str, data: bytes) -> str:
    if detect_domain is not None:
        try:
            return detect_domain(path, data)
        except Exception:
            pass
    m = detect_magic(data)
    if not m:
        return "unknown"
    if m.kind in IMAGE_KINDS:
        return "image"
    if m.kind in AUDIO_KINDS:
        return "audio"
    if m.kind in VIDEO_KINDS:
        return "video"
    if m.kind in {"pdf", "docx", "pptx", "xlsx", "epub"}:
        return "document"
    if m.kind in {"zip", "gzip", "xz", "bz2"}:
        return "archive"
    if m.kind == "text":
        return "text"
    return "other"

def preview_file(path: Optional[str], title: str = "Preview") -> None:
    st.markdown(f"#### {title}")
    if not path or not os.path.exists(path):
        st.info("No file available.")
        return
    size = os.path.getsize(path)
    data_head = Path(path).read_bytes()[:4096]
    m = detect_magic(Path(path).read_bytes()[: min(size, 1024 * 1024)]) if size <= 512 * 512 else detect_magic(data_head)
    ext = (m.ext if m else Path(path).suffix).lower()
    kind = m.kind if m else "unknown"
    #st.caption(f"{os.path.basename(path)} · {fmt_bytes(size)} · {kind}")

    try:
        if kind in IMAGE_KINDS or ext in {".png", ".jpg", ".jpeg", ".webp", ".gif", ".bmp", ".tif", ".tiff"}:
            st.image(
                path,
                use_container_width=IMAGE_PREVIEW_USE_CONTAINER_WIDTH,
                width=None if IMAGE_PREVIEW_USE_CONTAINER_WIDTH else IMAGE_PREVIEW_WIDTH,
            )
        elif kind in AUDIO_KINDS or ext in {".wav", ".mp3", ".ogg", ".opus", ".flac", ".m4a", ".aac"}:
            st.audio(path)
        elif kind in VIDEO_KINDS or ext in {".mp4", ".mov", ".mkv", ".avi", ".webm", ".m4v"}:
            st.video(path)
        elif kind == "text" or ext in {".txt", ".md", ".json", ".csv", ".tsv", ".log", ".xml", ".yaml", ".yml", ".html", ".py", ".js"}:
            text = Path(path).read_text(encoding="utf-8", errors="ignore")
            global _PREVIEW_FILE_CALL_COUNTER
            _PREVIEW_FILE_CALL_COUNTER += 1
            preview_key_src = f"{title}|{os.path.abspath(path)}|{size}|{_PREVIEW_FILE_CALL_COUNTER}"
            preview_key = "text_preview_" + hashlib.sha256(preview_key_src.encode("utf-8", errors="ignore")).hexdigest()[:16]
            st.text_area("Text preview", value=text[:15000], height=TEXT_PREVIEW_HEIGHT, label_visibility="collapsed", key=preview_key)
        elif kind == "pdf" or ext == ".pdf":
            st.info("PDF inline preview is not shown. Use download to inspect.")
        else:
            st.info("Preview is not available for this file type.")
    except Exception as e:
        st.warning(f"Preview failed: {e}")

def download_bytes_button(label: str, data: bytes, file_name: str, mime: str = "application/octet-stream") -> None:
    st.download_button(label, data=data, file_name=file_name, mime=mime, use_container_width=True)
