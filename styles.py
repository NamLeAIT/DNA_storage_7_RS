from __future__ import annotations

from ui_design_system.streamlit_style import apply_app_style


def apply_style() -> None:
    """Backward-compatible wrapper for older app.py imports."""
    apply_app_style()
