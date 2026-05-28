from __future__ import annotations

import streamlit as st

from config import APP_TITLE
from ui_design_system.streamlit_style import apply_app_style
from ui_design_system.ui_labels import PANEL_TITLES, STATUS
from panels import (
    render_panel_1_upload,
    render_panel_2_compression,
    render_panel_3_encoding,
    render_panel_4_experiment,
    render_panel_5_decoding,
    render_panel_6_analysis,
)


# The current UI has six visible panels. Error Simulation and Read Recovery are
# sections inside Strand Preparation, so the top stepper follows the rendered
# panel layout rather than the full conceptual pipeline.
APP_STEPS = [
    (1, PANEL_TITLES["input"]),
    (2, PANEL_TITLES["data_encoding"]),
    (3, PANEL_TITLES["dna_encoding"]),
    (4, PANEL_TITLES["strand_preparation"]),
    (5, PANEL_TITLES["file_decoding"]),
    (6, PANEL_TITLES["validation"]),
]


def _step_state(step_no: int) -> tuple[str, str]:
    checks = {
        1: bool(st.session_state.get("input_bytes")),
        2: bool(st.session_state.get("stored_bytes")),
        3: bool(st.session_state.get("dna")),
        4: bool(st.session_state.get("strand_rows")),
        5: bool(st.session_state.get("decoded_data")),
        6: bool(st.session_state.get("restored_info")),
    }
    if checks.get(step_no):
        return "done", STATUS["done"]
    previous_done = all(checks.get(i) for i in range(1, step_no)) if step_no > 1 else True
    if previous_done:
        return "current", "Next"
    return "", STATUS["waiting"]


def _render_hero() -> None:
    st.markdown(
        """
<div class="hero-card">
  <div class="hero-title">🧬 DNA Storage Pipeline</div>
  <div class="hero-subtitle">Compression-aware DNA encoding, sequencing simulation, read recovery, decoding, and validation.</div>
</div>
""",
        unsafe_allow_html=True,
    )


def _render_stepper() -> None:
    parts = ['<div class="pipeline-steps">']
    for n, label in APP_STEPS:
        css, state = _step_state(n)
        parts.append(
            f'<div class="pipeline-step {css}">'
            f'<div><span class="step-num">{n}</span><span class="step-name">{label}</span></div>'
            f'<div class="step-state">{state}</div>'
            f'</div>'
        )
    parts.append("</div>")
    st.markdown("".join(parts), unsafe_allow_html=True)


def render_app() -> None:
    st.set_page_config(page_title=APP_TITLE, page_icon="🧬", layout="wide")
    apply_app_style()
    _render_hero()
    _render_stepper()

    render_panel_1_upload()
    render_panel_2_compression()
    render_panel_3_encoding()
    render_panel_4_experiment()
    render_panel_5_decoding()
    render_panel_6_analysis()


if __name__ == "__main__":
    render_app()
