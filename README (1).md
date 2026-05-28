# Pipeline Structure

Use this structure for all future UI updates.

## 1. Input

Purpose: load the original file.

Show:
- input type
- detected type
- input size
- input preview
- input binary download

Use labels:
- `Input preview`
- `Input binary`
- `Download input binary`

---

## 2. Data Encoding

Purpose: decide which bytes will be stored in DNA.

Options:
- No compression
- Compression

For images with no compression:
- Original file bytes
- RGB pixels
- Grayscale pixels
- Binary image pixels

Main output:
- stored bytes
- stored binary bitstring

Use labels:
- `Storage method`
- `Stored size`
- `Stored type`
- `Run Data Encoding`
- `Download stored data`
- `Download stored binary`

---

## 3. DNA Encoding

Purpose: convert stored bytes into DNA.

Mapping methods:
- SM
- R∞
- R2
- R1
- R0
- Design Method A

Use labels:
- `DNA mapping`
- `DNA length`
- `GC content`
- `Longest HP`
- `Run DNA Encoding`
- `Download encoded DNA`
- `Download encoded binary`

---

## 4. Strand Preparation

Purpose: split encoded DNA into strand structures.

Structure:
```text
FBR + SI + Payload + Filler + RBR
```

Use labels:
- `Strand design`
- `Prepared strands`
- `Total strand length`
- `Inspect prepared strand`
- `Run Strand Preparation`
- `Download prepared strands`

Avoid:
- Clean strand as the main stage name
- Sample

---

## 5. Error Simulation

Purpose: add controlled errors after strand preparation and generate sequencing reads.

Subsections:
1. Strand-level Errors
2. Sequencing Read Errors

Use labels:
- `Strand-level errors`
- `Error target`
- `Substitution`
- `Insertion`
- `Deletion`
- `Run Add Errors`
- `Error strands`
- `Added errors`
- `Sequencing read errors`
- `Coverage`
- `Run Sequencing Simulation`
- `Sequencing reads`
- `Sequencing read errors`

---

## 6. Read Recovery

Purpose: recover DNA from noisy sequencing reads.

Use labels:
- `Read Recovery`
- `Run Read Recovery`
- `Recovered strands`
- `Reads recovered`
- `DNA ready for decoding`
- `Recovered DNA`
- `Download recovered DNA`
- `Download recovery table`

Meaning:
- `Reads recovered`: how well the sequencing reads were reconstructed back into the sequencing input.
- `DNA ready for decoding`: how closely recovered DNA matches the original encoded DNA.

---

## 7. File Decoding

Purpose: decode DNA into file bytes.

Use labels:
- `Input DNA`
- `Original encoded DNA`
- `Recovered DNA`
- `DNA mapping`
- `Run Decode`
- `Decoded size`
- `Restored type`
- `File can open`
- `Download decoded file`
- `Download decoded binary`

---

## 8. Validation

Purpose: check whether the decoded output is correct.

Use labels:
- `Restored type`
- `Restored size`
- `Restored correctly`
- `Restored preview`

For raw image pixels, also show:
- `Pixel accuracy`
- `Changed pixels`

Avoid:
- SHA-256 as a main UI metric
- File match
- Exact match
