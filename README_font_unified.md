# Font-unified app update

This package centralizes font family and font size into:

```text
ui_design_system/design_tokens.py
ui_design_system/streamlit_style.py
```

It uses a 3-size typography system:

```text
Large  = 20px
Medium = 15px
Small  = 12.5px
```

Replace these files in your app folder:

```text
app.py
panels.py
styles.py
ui_design_system/
```

`styles.py` is now only a compatibility wrapper. New code should import:

```python
from ui_design_system.streamlit_style import apply_app_style
```

The UI labels are still controlled by:

```text
ui_design_system/ui_labels.py
```

To rename `Design Method A`, edit:

```python
MAPPING_DISPLAY["New Design"] = "Protected Design"
```
