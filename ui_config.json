"""
Reusable Streamlit CSS for the DNA Storage app.

All text sizing must come from design_tokens.TYPOGRAPHY.
The app uses only three font sizes:
- Large: panel titles, hero title, metric values
- Medium: labels, buttons, normal UI text
- Small: captions, tables, DNA/binary previews
"""

import streamlit as st
from .design_tokens import COLORS, TYPOGRAPHY, RADIUS, SHADOWS, SPACING


def apply_app_style() -> None:
    """Apply global app styling."""
    lg = TYPOGRAPHY["large_size"]
    md = TYPOGRAPHY["medium_size"]
    sm = TYPOGRAPHY["small_size"]
    st.markdown(
        f"""
<style>
:root {{
    --bg: {COLORS["background"]};
    --surface: {COLORS["surface"]};
    --surface-soft: {COLORS["surface_soft"]};
    --text: {COLORS["text"]};
    --muted: {COLORS["muted"]};
    --border: {COLORS["border"]};
    --primary: {COLORS["primary"]};
    --primary-soft: {COLORS["primary_soft"]};
    --success-soft: {COLORS["success_soft"]};
    --warning-soft: {COLORS["warning_soft"]};
    --danger-soft: {COLORS["danger_soft"]};
    --font-main: {TYPOGRAPHY["font_family"]};
    --font-mono: {TYPOGRAPHY["mono_family"]};
    --text-lg: {lg};
    --text-md: {md};
    --text-sm: {sm};
}}

html, body, [class*="css"] {{
    font-family: var(--font-main);
    color: var(--text);
}}

.stApp {{
    background: linear-gradient(180deg, var(--bg) 0%, #EEF5FF 100%);
}}

.block-container {{
    padding-top: {SPACING["page_top"]};
    padding-bottom: {SPACING["page_bottom"]};
    max-width: 1280px;
}}

h1, h2 {{
    font-size: var(--text-lg) !important;
    font-weight: {TYPOGRAPHY["title_weight"]} !important;
    letter-spacing: -0.025em;
    color: var(--text);
}}

h3, h4 {{
    font-size: var(--text-md) !important;
    font-weight: {TYPOGRAPHY["section_title_weight"]} !important;
    color: var(--text);
}}

p, label, div, span {{
    font-size: var(--text-md);
}}

small, .small-note, .step-state, .step-num,
div[data-testid="stMetricLabel"] p,
div[data-testid="stDataFrame"] * {{
    font-size: var(--text-sm) !important;
}}

/* Main Streamlit bordered containers: soft cards */
div[data-testid="stVerticalBlockBorderWrapper"] {{
    background: linear-gradient(180deg, rgba(255,255,255,0.98), rgba(248,251,255,0.98));
    border: 1px solid var(--border) !important;
    border-radius: {RADIUS["card"]} !important;
    box-shadow: {SHADOWS["card"]};
}}

/* Metrics */
div[data-testid="stMetric"] {{
    background: rgba(255,255,255,0.90);
    border: 1px solid rgba(90,125,156,0.18);
    border-radius: {RADIUS["metric"]};
    padding: 0.48rem 0.65rem;
}}

div[data-testid="stMetricLabel"] p {{
    color: var(--muted) !important;
}}

div[data-testid="stMetricValue"] {{
    font-size: var(--text-lg) !important;
    font-weight: {TYPOGRAPHY["title_weight"]} !important;
    color: var(--text) !important;
}}

/* Buttons */
.stButton button, .stDownloadButton button {{
    border-radius: {RADIUS["button"]} !important;
    font-weight: {TYPOGRAPHY["strong_weight"]} !important;
    border: 1px solid rgba(37,99,235,0.20) !important;
    font-size: var(--text-md) !important;
}}

/* Inputs */
div[data-baseweb="input"] input,
div[data-baseweb="select"] div,
textarea {{
    border-radius: 11px !important;
    font-size: var(--text-md) !important;
}}

/* DNA/code previews */
code, pre, textarea {{
    font-family: var(--font-mono) !important;
}}

textarea, code, pre {{
    font-size: var(--text-sm) !important;
}}

/* Hero */
.hero-card {{
    padding: 1.05rem 1.15rem;
    border-radius: 20px;
    border: 1px solid rgba(37,99,235,0.16);
    background: radial-gradient(circle at 12% 10%, rgba(37,99,235,0.12), transparent 34%),
                linear-gradient(135deg, rgba(255,255,255,0.96), rgba(239,246,255,0.88));
    box-shadow: {SHADOWS["card"]};
    margin-bottom: 0.85rem;
}}
.hero-title {{
    font-size: var(--text-lg);
    font-weight: {TYPOGRAPHY["title_weight"]};
    letter-spacing: -0.025em;
    color: var(--text);
    margin-bottom: 0.22rem;
}}
.hero-subtitle {{
    color: var(--muted);
    font-size: var(--text-md);
    line-height: 1.45;
}}

/* Stepper */
.pipeline-steps {{
    display: grid;
    grid-template-columns: repeat(6, minmax(0, 1fr));
    gap: 0.5rem;
    margin: 0.8rem 0 1.0rem 0;
}}
.pipeline-step {{
    border: 1px solid rgba(90,125,156,0.18);
    border-radius: 14px;
    background: rgba(255,255,255,0.80);
    padding: 0.55rem 0.55rem;
    min-height: 58px;
}}
.pipeline-step.done {{
    background: linear-gradient(180deg, rgba(220,252,231,0.92), rgba(240,253,244,0.92));
    border-color: rgba(34,197,94,0.25);
}}
.pipeline-step.current {{
    background: linear-gradient(180deg, rgba(219,234,254,0.98), rgba(239,246,255,0.98));
    border-color: rgba(37,99,235,0.35);
}}
.step-num {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 20px;
    height: 20px;
    border-radius: 999px;
    background: rgba(37,99,235,0.10);
    color: var(--primary);
    font-weight: {TYPOGRAPHY["strong_weight"]};
    margin-right: 0.35rem;
}}
.step-name {{
    font-weight: {TYPOGRAPHY["strong_weight"]};
    color: var(--text);
    font-size: var(--text-sm);
}}
.step-state {{
    color: var(--muted);
    margin-top: 0.18rem;
}}

/* Section headers */
.step-heading {{
    display: flex;
    align-items: center;
    gap: 0.55rem;
    margin: 0.12rem 0 0.9rem 0;
}}
.step-badge {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 30px;
    height: 30px;
    border-radius: 10px;
    background: var(--primary-soft);
    color: var(--primary);
    font-weight: {TYPOGRAPHY["title_weight"]};
    font-size: var(--text-md);
}}
.step-title {{
    font-size: var(--text-lg);
    font-weight: {TYPOGRAPHY["title_weight"]};
    letter-spacing: -0.02em;
    color: var(--text);
}}

.pipeline-box {{
    padding: 0.75rem 0.9rem;
    border-radius: 15px;
    border: 1px solid rgba(90,125,156,0.18);
    background: rgba(255,255,255,0.76);
    color: var(--muted);
}}

/* DNA region visualization */
.region-tag {{
    display: inline-block;
    margin: 3px 4px 3px 0;
    padding: 4px 6px;
    border-radius: {RADIUS["tag"]};
    font-family: var(--font-mono);
    font-size: var(--text-sm);
    word-break: break-all;
}}
.error-base {{
    background: {COLORS["error_bg"]};
    color: {COLORS["error_text"]};
    font-weight: {TYPOGRAPHY["strong_weight"]};
    border-radius: 3px;
    padding: 0 2px;
}}

/* Tables */
div[data-testid="stDataFrame"] {{
    border-radius: 14px;
    overflow: hidden;
}}

@media (max-width: 900px) {{
    .pipeline-steps {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
}}
</style>
""",
        unsafe_allow_html=True,
    )
