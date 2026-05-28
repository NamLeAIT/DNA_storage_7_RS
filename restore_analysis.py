from __future__ import annotations

import gzip
import io
import lzma
import math
import os
import zipfile
from pathlib import Path
from typing import Any, Dict, Optional

try:
    from PIL import Image
except Exception:
    Image = None

try:
    import numpy as np
except Exception:
    np = None

from utils_bits_v2 import detect_magic, safe_basename, sha256_bytes
from ui_helpers import magic_dict
from config import IMAGE_KINDS


def extract_zip_first_file(data: bytes, out_dir: str) -> Optional[str]:
    try:
        with zipfile.ZipFile(io.BytesIO(data), "r") as zf:
            members = [m for m in zf.namelist() if not m.endswith("/")]
            if not members:
                return None
            member = members[0]
            safe = safe_basename(os.path.basename(member) or "extracted.bin")
            path = Path(out_dir) / safe
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(zf.read(member))
            return str(path)
    except Exception:
        return None

def _make_image_preview_png(path: str, out_dir: str, preferred_name: str) -> Optional[str]:
    """Create a browser/Streamlit-friendly PNG preview for any restored image."""
    if Image is None or not path or not os.path.exists(path):
        return None
    try:
        m = detect_magic(Path(path).read_bytes()[:1024 * 1024])
        if not m or m.kind not in IMAGE_KINDS:
            return None
        img = Image.open(path)
        # PNG previews are more reliable in Streamlit/browser than WEBP/BMP/TIFF.
        if img.mode not in {"RGB", "RGBA", "L"}:
            img = img.convert("RGBA" if "A" in img.getbands() else "RGB")
        preview = Path(out_dir) / f"{preferred_name}_preview.png"
        img.save(preview, format="PNG")
        return str(preview)
    except Exception:
        return None


def write_restored_file(data: bytes, out_dir: str, preferred_name: str = "restored") -> Dict[str, Any]:
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    m = detect_magic(data)
    ext = m.ext if m else ".bin"
    kind = m.kind if m else "unknown"
    file_path = Path(out_dir) / f"{preferred_name}{ext}"
    file_path.write_bytes(data)

    preview_path = str(file_path)
    extracted_path = None
    image_preview_path = None
    try:
        if kind == "zip" or kind in {"docx", "pptx", "xlsx", "epub"}:
            extracted_path = extract_zip_first_file(data, out_dir)
            preview_path = extracted_path or preview_path
        elif kind == "gzip":
            inner = gzip.decompress(data)
            inner_magic = detect_magic(inner)
            inner_ext = inner_magic.ext if inner_magic else ".bin"
            extracted_path = str(Path(out_dir) / f"{preferred_name}_gunzip{inner_ext}")
            Path(extracted_path).write_bytes(inner)
            preview_path = extracted_path
        elif kind == "xz":
            inner = lzma.decompress(data, format=lzma.FORMAT_XZ)
            inner_magic = detect_magic(inner)
            inner_ext = inner_magic.ext if inner_magic else ".bin"
            extracted_path = str(Path(out_dir) / f"{preferred_name}_unxz{inner_ext}")
            Path(extracted_path).write_bytes(inner)
            preview_path = extracted_path
        elif kind == "bz2":
            inner = bz2.decompress(data)
            inner_magic = detect_magic(inner)
            inner_ext = inner_magic.ext if inner_magic else ".bin"
            extracted_path = str(Path(out_dir) / f"{preferred_name}_bunzip2{inner_ext}")
            Path(extracted_path).write_bytes(inner)
            preview_path = extracted_path
    except Exception:
        # Keep native restored container if extraction fails.
        pass

    # Always create a PNG preview when the restored/extracted content is an image.
    # The stored byte stream is not changed; this is only for reliable UI display.
    image_preview_path = _make_image_preview_png(preview_path, out_dir, preferred_name)
    if image_preview_path:
        preview_path = image_preview_path

    return {
        "file_path": str(file_path),
        "preview_path": preview_path,
        "extracted_path": extracted_path,
        "image_preview_path": image_preview_path,
        "magic": magic_dict(data),
        "size_bytes": len(data),
        "sha256": sha256_bytes(data),
    }

def image_metrics(path_a: str, path_b: str) -> Dict[str, Any]:
    if Image is None or np is None:
        return {"ok": False, "reason": "PIL/numpy unavailable"}
    try:
        a = Image.open(path_a).convert("RGB")
        b = Image.open(path_b).convert("RGB")
        if a.size != b.size:
            b = b.resize(a.size)
        arr_a = np.asarray(a).astype("float32")
        arr_b = np.asarray(b).astype("float32")
        mse = float(np.mean((arr_a - arr_b) ** 2))
        mae = float(np.mean(np.abs(arr_a - arr_b)))
        psnr = 99.0 if mse <= 1e-12 else float(20.0 * math.log10(255.0 / math.sqrt(mse)))

        # Try scikit-image SSIM; otherwise use a simple global approximation.
        ssim_val = None
        try:
            from skimage.metrics import structural_similarity as ssim
            ssim_val = float(ssim(arr_a.astype("uint8"), arr_b.astype("uint8"), channel_axis=2, data_range=255))
        except Exception:
            x = arr_a.reshape(-1, 3)
            y = arr_b.reshape(-1, 3)
            C1 = (0.01 * 255) ** 2
            C2 = (0.03 * 255) ** 2
            vals = []
            for ch in range(3):
                xx = x[:, ch]
                yy = y[:, ch]
                mux, muy = float(xx.mean()), float(yy.mean())
                vx, vy = float(xx.var()), float(yy.var())
                cov = float(((xx - mux) * (yy - muy)).mean())
                vals.append(((2 * mux * muy + C1) * (2 * cov + C2)) / ((mux * mux + muy * muy + C1) * (vx + vy + C2)))
            ssim_val = float(np.mean(vals))

        return {"Validation": True, "mse": mse, "mae": mae, "psnr": psnr, "ssim": ssim_val, "size": a.size}
    except Exception as e:
        return {"Validation": False, "reason": str(e)}

def text_similarity(path_a: str, path_b: str, cap: int = 20000) -> Dict[str, Any]:
    try:
        a = Path(path_a).read_text(encoding="utf-8", errors="ignore")[:cap]
        b = Path(path_b).read_text(encoding="utf-8", errors="ignore")[:cap]
        exact = a == b
        # Simple character accuracy by positional comparison, not full Levenshtein.
        L = max(1, max(len(a), len(b)))
        same = sum(1 for x, y in zip(a, b) if x == y)
        return {"Validation": True, "exact": exact, "char_position_accuracy": same / L, "len_a": len(a), "len_b": len(b)}
    except Exception as e:
        return {"Validation": False, "reason": str(e)}
