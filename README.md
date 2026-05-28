# Compression-aware DNA Storage App

Pipeline:
1. Upload
2. Encode — Bytes (No compression or Compression)
3. Encode — DNA (SM, R∞, R2, R1, R0, Design Method A)
4. DNA Strand Prep + Advanced Add Errors + Wet-lab Read Simulation + Read Reconstruction
5. Decode
6. Validate

Key UI changes:
- New Design is displayed as Design Method A.
- All processing buttons use the Run prefix.
- Wet-lab controls are simplified.
- Each stage provides downloads for the relevant bytes, binary, DNA, strands, reads, errors, or decoded output.
