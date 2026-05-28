from __future__ import annotations

import tempfile
from pathlib import Path

APP_TITLE = "Compression-Aware DNA Storage Pipeline"

# All temporary working files are written here.
WORK_ROOT = Path(tempfile.gettempdir()) / "headerless_dna_pipeline"
WORK_ROOT.mkdir(parents=True, exist_ok=True)

# Preview / layout controls
IMAGE_PREVIEW_USE_CONTAINER_WIDTH = False
IMAGE_PREVIEW_WIDTH = 200
TEXT_PREVIEW_HEIGHT = 250
DNA_PREVIEW_HEIGHT = 150
FRAGMENT_PREVIEW_HEIGHT = 180

SELF_DESCRIBING_KINDS = {
    "png", "jpeg", "webp", "gif", "bmp", "tiff",
    "pdf",
    "zip", "docx", "pptx", "xlsx", "epub",
    "gzip", "xz", "bz2",
    "mp4", "avi", "mkv_webm",
    "wav", "mp3", "flac", "ogg", "opus_ogg",
    "text",
}

IMAGE_KINDS = {"png", "jpeg", "webp", "gif", "bmp", "tiff"}
AUDIO_KINDS = {"wav", "mp3", "flac", "ogg", "opus_ogg"}
VIDEO_KINDS = {"mp4", "avi", "mkv_webm"}

MAPPING_OPTIONS = [
    "Simple Mapping",
    "RINF_B16",
    "R2_B15",
    "R1_B12",
    "R0_B9",
    "New Design",
    "Reed-Solomon",
]
