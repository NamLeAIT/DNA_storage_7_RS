"""
Design tokens for the DNA Storage Streamlit app.

Edit colors, font sizes, spacing, and radius here.
Use these constants in style/helper modules instead of hard-coding values.
"""

# -----------------------------------------------------------------------------
# Color palette
# -----------------------------------------------------------------------------

COLORS = {
    # Core
    "background": "#F8FAFC",
    "surface": "#FFFFFF",
    "surface_soft": "#F1F5F9",
    "text": "#0F172A",
    "muted": "#64748B",
    "border": "#D8E1EC",

    # Brand
    "primary": "#2563EB",
    "primary_soft": "#DBEAFE",
    "primary_text": "#1E3A8A",

    # Status
    "success": "#16A34A",
    "success_soft": "#DCFCE7",
    "warning": "#D97706",
    "warning_soft": "#FEF3C7",
    "danger": "#DC2626",
    "danger_soft": "#FEE2E2",

    # DNA region colors
    "fbr_bg": "#DBEAFE",
    "fbr_text": "#1E3A8A",
    "si_bg": "#EDE9FE",
    "si_text": "#4C1D95",
    "payload_bg": "#DCFCE7",
    "payload_text": "#14532D",
    "filler_bg": "#F1F5F9",
    "filler_text": "#475569",
    "rbr_bg": "#FFEDD5",
    "rbr_text": "#7C2D12",
    "error_bg": "#FECACA",
    "error_text": "#7F1D1D",
}

REGION_COLORS = {
    "FBR": (COLORS["fbr_bg"], COLORS["fbr_text"]),
    "SI": (COLORS["si_bg"], COLORS["si_text"]),
    "Payload": (COLORS["payload_bg"], COLORS["payload_text"]),
    "Filler": (COLORS["filler_bg"], COLORS["filler_text"]),
    "RBR": (COLORS["rbr_bg"], COLORS["rbr_text"]),
    "Error": (COLORS["error_bg"], COLORS["error_text"]),
}

# -----------------------------------------------------------------------------
# Typography
# -----------------------------------------------------------------------------

TYPOGRAPHY = {
    "font_family": 'Inter, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
    "mono_family": '"SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace',

    # Three-size typography system only.
    # Large  = panel titles, hero title, metric values
    # Medium = labels, buttons, body text, section names
    # Small  = captions, tables, DNA/binary previews
    "large_size": "20px",
    "medium_size": "15px",
    "small_size": "15px",

    # Backward-compatible aliases used by existing CSS helpers.
    "title_size": "20px",
    "panel_title_size": "20px",
    "section_title_size": "15px",
    "body_size": "15px",
    "metric_size": "20px",

    "title_weight": 760,
    "panel_title_weight": 760,
    "section_title_weight": 700,
    "body_weight": 400,
    "medium_weight": 560,
    "strong_weight": 700,
}

# -----------------------------------------------------------------------------
# Layout
# -----------------------------------------------------------------------------

SPACING = {
    "page_top": "1.0rem",
    "page_bottom": "2.4rem",
    "card_padding": "1.0rem",
    "small_padding": "0.5rem",
    "gap": "1rem",
}

RADIUS = {
    "card": "15px",
    "button": "15px",
    "metric": "15px",
    "tag": "15px",
}

SHADOWS = {
    "card": "0 10px 30px rgba(15, 23, 42, 0.045)",
}

# -----------------------------------------------------------------------------
# Display limits
# -----------------------------------------------------------------------------

DISPLAY = {
    "dna_preview_chars": 600,
    "binary_preview_chars": 3000,
    "sequence_preview_chars": 900,
    "image_preview_width": 220,
    "table_height": 320,
}
