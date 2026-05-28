# DNA Storage App UI Design System

This folder defines the visual style, naming rules, UI labels, panel structure, and reusable helpers for the DNA Storage Streamlit app.

Goal:

```text
One pipeline, one naming system, one visual style.
```

Use this package when you want to update the app UI without editing labels, colors, and layout rules manually across many files.

---

## Recommended app flow

```text
Input
→ Data Encoding
→ DNA Encoding
→ Strand Preparation
→ Error Simulation
→ Read Recovery
→ File Decoding
→ Validation
```

This naming separates the pipeline clearly:

| Stage | Meaning |
|---|---|
| Input | Upload file and inspect original data |
| Data Encoding | Choose no-compression/compression and create stored bytes |
| DNA Encoding | Convert stored bytes to DNA using SM/R∞/R2/R1/R0/Design Method A |
| Strand Preparation | Split DNA into strands with FBR/SI/Payload/Filler/RBR |
| Error Simulation | Add visible strand errors and sequencing read errors |
| Read Recovery | Reconstruct DNA from noisy sequencing reads |
| File Decoding | Decode recovered DNA back into bytes/file |
| Validation | Check whether restored data is correct |

---

## How to use in the app

Suggested imports in `panels.py`:

```python
from ui_design_system.ui_labels import LABELS, PANEL_TITLES, BUTTONS, METRICS
from ui_design_system.design_tokens import COLORS, TYPOGRAPHY, SPACING
from ui_design_system.streamlit_style import apply_app_style
```

At app startup:

```python
apply_app_style()
```

Then replace hard-coded strings like:

```python
step_header(4, "DNA Strand Preparation")
```

with:

```python
step_header(4, PANEL_TITLES["strand_preparation"])
```

---

## Main rule

Avoid scattered UI strings in `panels.py`.

Good:

```python
st.button(BUTTONS["run_decode"])
```

Bad:

```python
st.button("Decoding")
```

