# Font Guidelines

The app uses only three text sizes.

| Level | Token | Default | Used for |
|---|---|---:|---|
| Large | `large_size` | 20px | panel titles, hero title, metric values |
| Medium | `medium_size` | 15px | labels, buttons, normal text, section names |
| Small | `small_size` | 12.5px | captions, tables, DNA/binary previews |

Rules:

1. Do not hard-code `font-size` in `app.py` or `panels.py`.
2. Do not define new font families outside `design_tokens.py`.
3. Use the normal UI font for dashboard text.
4. Use the mono font only for DNA, binary, code, and sequence previews.
5. If font size needs to change, edit `design_tokens.py` only.

Main files:

```text
ui_design_system/design_tokens.py      # font sizes, colors, spacing
ui_design_system/streamlit_style.py    # CSS that uses the tokens
app.py                                 # no font styling
panels.py                              # no font-size styling
```
