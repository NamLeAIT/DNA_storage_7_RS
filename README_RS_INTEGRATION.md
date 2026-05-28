# Reed-Solomon integration update

This version updates the main app with the new Reed-Solomon design and a vertical Panel 4 flow.

## Main updates

- Added `Reed-Solomon` to DNA Encoding.
- Replaced the previous RS/demo-style implementation with a dna_rs_coding-style concatenated Reed-Solomon design:
  - Inner RS: GF(2^6), 1 symbol = 6 bits = 3 nt.
  - Outer RS: GF(2^14), across DNA strands.
- Reed-Solomon strand rows are used directly in Strand Preparation.
- Panel 4 no longer uses tabs for wet-lab/read recovery.
- Panel 4 now flows vertically:
  1. Strand Preparation
  2. Optional strand-level errors
  3. Sequencing Simulation
  4. Read Recovery
- Strand-level error settings are hidden behind a checkbox.
- Sequencing Simulation settings remain available.
- Read Recovery runs after Sequencing Simulation and shows results below.

## Notes

- Substitution errors are best supported.
- Light insertion/deletion can be handled through the added alignment step before RS correction.
- Heavy insertion/deletion may still cause partial recovery or failure.
