# Final RS Strand Design Update

This app version applies the final requested design updates.

## Main changes

1. Mapping list keeps:
   - SM
   - R∞
   - R2
   - R1
   - R0
   - Protected Design
   - Reed-Solomon

2. DNA Decoding is symmetric with DNA Encoding:
   - The selected DNA mapping and saved codec metadata are reused for decoding.

3. Reed-Solomon strand structure is now full 125 nt:
   - FBR + Index + Payload + RS parity + RBR
   - FBR = 20 nt
   - Index = 12 nt
   - Payload = 48 nt
   - RS parity = 24 nt
   - RBR = 21 nt
   - Total = 125 nt

4. Reed-Solomon index is masked/whitened:
   - Avoids AAAAAAAAAAAA for index 0.
   - Reduces bad homopolymer/GC behavior.
   - Index remains reversible for decoding.

5. Normal method strand structure stays:
   - FBR + SI + Payload + Filler + RBR

6. UI wording:
   - Uses Strand Design rather than Strand Preparation.

7. Fair expansion metrics:
   - DNA expansion in DNA Encoding compares coding-level DNA only.
   - Strand Design expansion compares full sequencing-ready strands.

## Test summary

The following core tests passed:
- Reed-Solomon core DNA decode.
- Reed-Solomon full 125-nt strand decode.
- Sequencing simulation + read recovery with 1% substitution.
- Light insertion/deletion simulation + read recovery.
- Protected Design basic encode/decode smoke test.
