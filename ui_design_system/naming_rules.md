# Naming Rules

## General

Use short labels. Avoid long explanatory subtitles in the main UI.

Good:
```text
Run Decode
Restored correctly
DNA ready for decoding
```

Avoid:
```text
The reconstructed DNA is not identical to the original encoded DNA...
Decode DNA acc.
SHA-256
```

---

## Buttons

All process buttons must start with `Run`.

Good:
```text
Run Data Encoding
Run DNA Encoding
Run Strand Preparation
Run Add Errors
Run Sequencing Simulation
Run Read Recovery
Run Decode
```

Avoid:
```text
Decoding
Start Encoding
Generate noisy reads and reconstruct
Run Simulation
```

---

## DNA terminology

Use:
```text
DNA
base string
prepared strand
error strand
sequencing reads
read errors
recovered DNA
```

Avoid:
```text
sequencings
sample
clean strand as a main panel name
```

`Clean strand` can be used only when comparing clean vs error strand visually.

---

## Mapping names

Use public names:

```text
SM
R∞
R2
R1
R0
Design Method A
```

Internal key can remain:

```text
New Design
```

but UI must display:

```text
Design Method A
```

---

## Recovery terms

Use:
```text
Read Recovery
Reads recovered
DNA ready for decoding
Restored correctly
```

Avoid:
```text
ECC reconstruction
Decode DNA acc.
Payload accuracy
File match
```

Explanation:
- Read Recovery is consensus/read-level correction.
- Design Method A performs local block-level repair during decoding.
- Validation checks final data correctness.
