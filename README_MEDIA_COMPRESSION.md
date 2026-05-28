# Media compression update

This version enables audio/video compression candidates in Streamlit.

## Added audio candidates

- FLAC lossless
- WAV PCM
- OGG/Opus
- OGG/Vorbis
- MP3
- AAC/M4A

## Added video candidates

- H.264 MP4
- H.265 MP4 when available
- VP9 WebM
- AV1 WebM when available
- H.264 MKV

## Streamlit Cloud

The app includes:

- `requirements.txt` with `imageio-ffmpeg`
- `packages.txt` with `ffmpeg`

`compression_pipeline.py` now uses FFmpeg from PATH first, then falls back to `imageio_ffmpeg.get_ffmpeg_exe()`.
