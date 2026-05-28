"""
Example usage in panels.py.

This file is not required by the app. It shows how to replace hard-coded UI labels.
"""

import streamlit as st
from ui_design_system.ui_labels import PANEL_TITLES, BUTTONS, METRICS, display_mapping


def example_decode_panel(mapping: str, valid: bool) -> None:
    st.header(PANEL_TITLES["file_decoding"])

    c1, c2, c3 = st.columns(3)
    c1.metric(METRICS["dna_mapping"], display_mapping(mapping))
    c2.metric(METRICS["file_can_open"], "Yes" if valid else "No")
    c3.metric(METRICS["restored_correctly"], "Yes")

    if st.button(BUTTONS["run_decode"]):
        st.write("Run decoder here.")
