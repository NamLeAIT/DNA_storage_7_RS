# Migration Checklist for panels.py

Use this checklist when updating the current app.

## Replace panel names

- `Upload` → `Input`
- `Encode — Bytes` → `Data Encoding`
- `Encode — DNA` → `DNA Encoding`
- `DNA Strand Preparation` → `Strand Preparation`
- `Advanced Setting` → `Strand-level Errors`
- `Wet-lab Sequencing Simulation` → `Sequencing Read Errors`
- `Recovery result from error strand` → `Read Recovery`
- `Decode` → `File Decoding`
- `Validate` → `Validation`

## Replace buttons

- `Run Compression` → `Run Data Encoding`
- `Run Encode` → `Run DNA Encoding`
- `Run DNA Strand Preparation` → `Run Strand Preparation`
- `Run Adding Errors` → `Run Add Errors`
- `Run Simulation` → `Run Sequencing Simulation`
- `Run Strand Recovery` → `Run Read Recovery`
- `Decoding` → `Run Decode`

## Replace metrics

- `Clean strands` → `Prepared strands`
- `Total strand nt` → `Total strand length`
- `Number of error strands` → `Error strands`
- `Number of added errors` → `Added errors`
- `Read errors` → `Sequencing read errors`
- `DNA ready for decode` → `DNA ready for decoding`
- `File extension` → `Restored type`
- `File recovery` → `File can open`
- `Exact match` / `File match` → `Restored correctly`

## Download labels

- `Download compressed data` → `Download stored data`
- `Download binary string` → `Download stored binary`
- `Download clean strands` → `Download prepared strands`
- `Download reconstruction table` → `Download recovery table`

## Delete or hide

- `_toolkit_rows_from_meta()` if Toolkit RS is not part of the active pipeline.
- Long captions under main controls.
- SHA-256 from the main UI.
- Read reconstruction method selector if the default is always index-aware.
- Seed from the main UI unless debugging.
