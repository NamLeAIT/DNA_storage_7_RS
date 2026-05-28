"""
Central UI labels for the DNA Storage app.

Edit this file when you want to rename panels, buttons, metrics, tabs, or mapping names.
Do not hard-code repeated UI strings in panels.py.
"""

# -----------------------------------------------------------------------------
# Panel titles
# -----------------------------------------------------------------------------

PANEL_TITLES = {
    "input": "Input",
    "data_encoding": "Compression",
    "dna_encoding": "DNA Encoding",
    "strand_preparation": "Strand Design",
    "error_simulation": "Error Simulation",
    "read_recovery": "Read Recovery",
    "file_decoding": "File Decoding",
    "validation": "Validation",
}

STEPPER = [
    PANEL_TITLES["input"],
    PANEL_TITLES["data_encoding"],
    PANEL_TITLES["dna_encoding"],
    PANEL_TITLES["strand_preparation"],
    PANEL_TITLES["error_simulation"],
    PANEL_TITLES["read_recovery"],
    PANEL_TITLES["file_decoding"],
    PANEL_TITLES["validation"],
]

# -----------------------------------------------------------------------------
# Buttons
# -----------------------------------------------------------------------------

BUTTONS = {
    "run_data_encoding": "Run Compression",
    "run_dna_encoding": "Run DNA Encoding",
    "run_strand_preparation": "Run Strand Design",
    "run_add_errors": "Run Add Errors",
    "run_sequencing_simulation": "Run Sequencing Simulation",
    "run_read_recovery": "Run Read Recovery",
    "run_decode": "Run Decode",
    "run_validation": "Run Validation",

    "download_input_binary": "Download input binary",
    "download_stored_data": "Download stored data",
    "download_stored_binary": "Download stored binary",
    "download_encoded_dna": "Download encoded DNA",
    "download_encoded_binary": "Download encoded binary",
    "download_prepared_strands": "Download designed strands",
    "download_error_strands": "Download error strands",
    "download_error_table": "Download error table",
    "download_noisy_reads": "Download noisy reads",
    "download_read_errors": "Download read errors",
    "download_recovered_dna": "Download recovered DNA",
    "download_recovery_table": "Download recovery table",
    "download_decoded_file": "Download decoded file",
    "download_decoded_binary": "Download decoded binary",
}

# -----------------------------------------------------------------------------
# Tabs
# -----------------------------------------------------------------------------

TABS = {
    "input_strands": "Input strands",
    "read_errors": "Read errors",
    "recovery_result": "Recovery result",
    "strand_level_errors": "Strand-level errors",
    "sequencing_errors": "Sequencing errors",
}

# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------

METRICS = {
    "input_type": "Input type",
    "detected_type": "Detected type",
    "input_size": "Input size",

    "storage_method": "Storage method",
    "stored_size": "Stored size",
    "stored_type": "Stored type",

    "dna_mapping": "DNA mapping",
    "dna_length": "DNA length",
    "dna_expansion": "DNA length increase",
    "gc_content": "GC content",
    "longest_hp": "Longest HP",

    "prepared_strands": "Designed strands",
    "total_strand_length": "Total strand length",
    "strand_design_expansion": "Strand Design length increase",

    "error_strands": "Error strands",
    "added_errors": "Added errors",

    "sequencing_reads": "Sequencing reads",
    "sequencing_read_errors": "Sequencing read errors",

    "recovered_strands": "Recovered strands",
    "reads_recovered": "Reads recovered",
    "dna_ready_for_decoding": "DNA ready for decoding",

    "input_dna": "Input DNA",
    "decoded_size": "Decoded size",
    "restored_type": "Restored type",
    "file_can_open": "File can open",
    "restored_correctly": "Restored correctly",

    "pixel_accuracy": "Pixel accuracy",
    "changed_pixels": "Changed pixels",
}

# -----------------------------------------------------------------------------
# Mapping display names
# -----------------------------------------------------------------------------

MAPPING_DISPLAY = {
    "Simple Mapping": "SM",
    "RINF_B16": "R∞",
    "R2_B15": "R2",
    "R1_B12": "R1",
    "R0_B9": "R0",
    "New Design": "Protected Design",
    "Reed-Solomon": "Reed-Solomon",
}

MAPPING_ORDER = [
    "Simple Mapping",
    "RINF_B16",
    "R2_B15",
    "R1_B12",
    "R0_B9",
    "New Design",
    "Reed-Solomon",
]

# -----------------------------------------------------------------------------
# Data source names
# -----------------------------------------------------------------------------

DATA_SOURCES = {
    "no_compression": "No compression",
    "compression": "Compression",
    "original_file_bytes": "Original file bytes",
    "rgb_pixels": "RGB pixels",
    "grayscale_pixels": "Grayscale pixels",
    "binary_image_pixels": "Binary image pixels",
}

# -----------------------------------------------------------------------------
# DNA region names
# -----------------------------------------------------------------------------

DNA_REGIONS = {
    "fbr": "FBR",
    "si": "SI",
    "payload": "Payload",
    "filler": "Filler",
    "rbr": "RBR",
}

# -----------------------------------------------------------------------------
# Status text
# -----------------------------------------------------------------------------

STATUS = {
    "yes": "Yes",
    "no": "No",
    "success": "Success",
    "fail": "Fail",
    "waiting": "Waiting",
    "ready": "Ready",
    "done": "Done",
}

def display_mapping(mapping: str) -> str:
    """Return user-facing name for an internal mapping key."""
    return MAPPING_DISPLAY.get(str(mapping), str(mapping))


# -----------------------------------------------------------------------------
# Field labels, messages, and download filenames
# -----------------------------------------------------------------------------

FIELDS = {
    "input_preview": "Input preview",
    "input_binary": "Input binary",
    "binary_bitstream": "Binary bitstream",
    "storage_method": "Storage method",
    "image_data_source": "Image data source",
    "binary_threshold": "Binary threshold",
    "pixel_format": "Pixel format",
    "image_size": "Image size",
    "raw_bytes": "Raw bytes",
    "pixel_view_used": "Pixel view used for DNA encoding",
    "compressed_output": "Compressed output",
    "dna_mapping": "DNA mapping",
    "base_string": "Base string",
    "strand_design": "Strand design",
    "total_strand_length": "Total strand length",
    "strand_design_expansion": "Strand Design length increase",
    "si_length": "SI length",
    "fbr": "FBR",
    "rbr": "RBR",
    "inspect_prepared_strand": "Inspect designed strand",
    "prepared_strand": "Designed strand",
    "clean_strand": "Clean strand",
    "strand_level_errors": "Strand-level Errors",
    "error_settings": "Error settings",
    "error_target": "Error target",
    "substitution": "Substitution Error",
    "insertion": "Insertion Error",
    "deletion": "Deletion Error",
    "seed": "Seed",
    "allow_indels": "Allow indels",
    "inspect_error_strand": "Inspect error strand",
    "error_strand": "Error strand",
    "sequencing_read_errors": "Sequencing Read Errors",
    "sequencing_error_settings": "Sequencing error settings",
    "coverage": "Reads per strand",
    "dropout": "Lost strands",
    "sequencing_input": "Sequencing input",
    "recovered_dna": "Recovered DNA",
    "input_dna_preview": "Input DNA preview",
    "restored_preview": "Restored preview",
    "image_comparison": "Image comparison",
    "text_comparison": "Text comparison",
}

MESSAGES = {
    "upload_to_start": "Upload a file to start.",
    "upload_first": "Upload a file first.",
    "choose_storage": "Choose storage method and run Data Encoding.",
    "run_data_encoding_first": "Run Data Encoding first.",
    "run_dna_encoding_first": "Run DNA Encoding first.",
    "run_strand_preparation": "Run Strand Design to continue.",
    "run_sequencing_simulation": "Run Sequencing Simulation.",
    "run_sequencing_first": "Run Sequencing Simulation first.",
    "run_read_recovery": "Run Read Recovery.",
    "run_decode_first": "Run Decode first.",
}

DOWNLOAD_FILES = {
    "input_binary": "input_binary.txt",
    "pixel_view_png": "input_pixel_view.png",
    "stored_data": "stored_data.bin",
    "stored_binary": "stored_binary.txt",
    "encoded_dna": "encoded_dna.txt",
    "encoded_binary": "encoded_binary.txt",
    "prepared_strands": "designed_strands.csv",
    "error_strands": "error_strands.csv",
    "error_table": "strand_error_table.csv",
    "sequencing_input_strands": "sequencing_input_strands.csv",
    "noisy_reads": "noisy_reads.csv",
    "read_errors": "sequencing_read_errors.csv",
    "recovered_dna": "recovered_dna.txt",
    "recovery_table": "recovery_table.csv",
    "decoded_binary": "decoded_binary.txt",
    "decoded_raw_pixels": "decoded_raw_pixels.bin",
    "decoded_raw_pixel_binary": "decoded_raw_pixel_binary.txt",
}

# Add any metric keys that older UI files may not include.
METRICS.update({
    "input_type": "Input type",
    "detected_type": "Detected type",
    "input_size": "Input size",
    "storage_method": "Storage method",
    "stored_size": "Stored size",
    "stored_type": "Stored type",
    "dna_mapping": "DNA mapping",
    "gc_content": "GC content",
    "prepared_strands": "Designed strands",
    "total_strand_length": "Total strand length",
    "strand_design_expansion": "Strand Design length increase",
    "error_strands": "Error strands",
    "added_errors": "Added errors",
    "sequencing_reads": "Sequencing reads",
    "sequencing_read_errors": "Sequencing read errors",
    "recovered_strands": "Recovered strands",
    "input_dna": "Input DNA",
    "restored_size": "Restored size",
})


# Reed-Solomon / vertical experiment additions
FIELDS.update({
    "add_strand_errors": "Advanced setting",
    "rs_design": "Reed-Solomon design",
    "sequencing_result": "Sequencing result",
    "recovery_result": "Recovery result",
})


# Final RS / Strand Design metric labels
METRICS.update({
    "dna_expansion": "DNA length increase",
    "strand_design_expansion": "Strand Design length increase",
    "prepared_strands": "Designed strands",
    "total_strand_length": "Total strand length",
})
FIELDS.update({
    "add_strand_errors": "Advanced setting",
    "sequencing_result": "Sequencing result",
    "recovery_result": "Recovery result",
})
