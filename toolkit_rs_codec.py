from __future__ import annotations

"""
Built-in DNAStorageToolkit-style baseline codec.

This module mirrors the structural codec used in DNAStorageToolkit's
1-encoding-decoding baseline:

    bytes -> 16-bit symbols -> 14 Reed-Solomon rows -> 120-nt strands

Each strand is:

    8 nt column index + 14 × 8 nt row symbols = 120 nt

The official toolkit uses Schifra Reed-Solomon over 16-bit symbols and a
fixed code length of 2^16-1. For a dependency-light Streamlit app, this module
implements a compatible systematic RS-style MDS erasure code over GF(2^16)
using the same 16-bit symbol and 2-bit DNA representation. It is intended as a
research baseline and app-side codec; the app also keeps an external-repo hook
for running the official C++ toolkit when available.
"""

import math
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

BASES = "ACGT"  # official Toolkit order: A=0, C=1, G=2, T=3
BASE2VAL = {b: i for i, b in enumerate(BASES)}
RS_ROWS = 14
INDEX_NT = 8
SYMBOL_NT = 8
STRAND_NT = INDEX_NT + RS_ROWS * SYMBOL_NT  # 120 nt
GF_SIZE = 1 << 16
GF_MASK = GF_SIZE - 1
PRIMITIVE_POLY = 0x1100B  # primitive for GF(2^16); alpha=2 has period 65535

_EXP: Optional[List[int]] = None
_LOG: Optional[List[int]] = None


def clean_dna(seq: str) -> str:
    return "".join(ch for ch in str(seq or "").upper() if ch in BASE2VAL)


def _init_tables() -> Tuple[List[int], List[int]]:
    global _EXP, _LOG
    if _EXP is not None and _LOG is not None:
        return _EXP, _LOG
    exp = [0] * (2 * (GF_SIZE - 1))
    log = [0] * GF_SIZE
    x = 1
    for i in range(GF_SIZE - 1):
        exp[i] = x
        log[x] = i
        x <<= 1
        if x & GF_SIZE:
            x ^= PRIMITIVE_POLY
        x &= GF_MASK
    for i in range(GF_SIZE - 1, 2 * (GF_SIZE - 1)):
        exp[i] = exp[i - (GF_SIZE - 1)]
    _EXP, _LOG = exp, log
    return exp, log


def gf_add(a: int, b: int) -> int:
    return (int(a) ^ int(b)) & GF_MASK


def gf_sub(a: int, b: int) -> int:
    return gf_add(a, b)  # characteristic two


def gf_mul(a: int, b: int) -> int:
    a &= GF_MASK
    b &= GF_MASK
    if a == 0 or b == 0:
        return 0
    exp, log = _init_tables()
    return exp[log[a] + log[b]]


def gf_div(a: int, b: int) -> int:
    a &= GF_MASK
    b &= GF_MASK
    if b == 0:
        raise ZeroDivisionError("GF division by zero")
    if a == 0:
        return 0
    exp, log = _init_tables()
    return exp[(log[a] - log[b]) % (GF_SIZE - 1)]


def gf_inv(a: int) -> int:
    a &= GF_MASK
    if a == 0:
        raise ZeroDivisionError("GF inverse of zero")
    exp, log = _init_tables()
    return exp[(GF_SIZE - 1 - log[a]) % (GF_SIZE - 1)]


def gf_pow_alpha(i: int) -> int:
    exp, _log = _init_tables()
    return exp[int(i) % (GF_SIZE - 1)]


def _symbol_to_dna(symbol: int) -> str:
    symbol &= GF_MASK
    # official order: high byte then low byte, MSB-first in 2-bit chunks
    out = []
    for shift in range(14, -1, -2):
        out.append(BASES[(symbol >> shift) & 0b11])
    return "".join(out)


def _dna_to_symbol(seq8: str) -> int:
    seq8 = clean_dna(seq8)
    if len(seq8) != 8:
        raise ValueError("A 16-bit Toolkit symbol must contain exactly 8 DNA bases")
    x = 0
    for ch in seq8:
        x = (x << 2) | BASE2VAL[ch]
    return x & GF_MASK


def _index_to_dna(index: int) -> str:
    if not (0 <= int(index) <= GF_MASK):
        raise ValueError("Toolkit column index must fit in 16 bits")
    return _symbol_to_dna(int(index))


def _dna_to_index(seq8: str) -> int:
    return _dna_to_symbol(seq8)


def _bytes_to_symbols(data: bytes) -> List[int]:
    out: List[int] = []
    for i in range(0, len(data), 2):
        hi = data[i]
        lo = data[i + 1] if i + 1 < len(data) else 0
        out.append(((hi << 8) | lo) & GF_MASK)
    return out


def _symbols_to_bytes(symbols: Sequence[int], file_size: int) -> bytes:
    out = bytearray()
    for s in symbols:
        s = int(s) & GF_MASK
        out.append((s >> 8) & 0xFF)
        out.append(s & 0xFF)
    if file_size >= 0:
        return bytes(out[: int(file_size)])
    return bytes(out)


def _column_points(n: int) -> List[int]:
    if n >= GF_SIZE:
        raise ValueError("Toolkit-style GF(2^16) codec supports at most 65,535 columns")
    return [gf_pow_alpha(i) for i in range(n)]


def _barycentric_weights(xs: Sequence[int]) -> List[int]:
    n = len(xs)
    weights: List[int] = []
    for i in range(n):
        denom = 1
        xi = xs[i]
        for j in range(n):
            if i == j:
                continue
            denom = gf_mul(denom, gf_sub(xi, xs[j]))
        weights.append(gf_inv(denom))
    return weights


def _lagrange_eval(xs: Sequence[int], ys: Sequence[int], x: int, weights: Optional[Sequence[int]] = None) -> int:
    # Return exact known value if x is one of the sample points.
    for xi, yi in zip(xs, ys):
        if int(xi) == int(x):
            return int(yi) & GF_MASK
    if weights is None:
        weights = _barycentric_weights(xs)
    prod = 1
    for xi in xs:
        prod = gf_mul(prod, gf_sub(x, xi))
    total = 0
    for xi, yi, wi in zip(xs, ys, weights):
        term = gf_mul(int(yi) & GF_MASK, wi)
        term = gf_div(term, gf_sub(x, xi))
        total ^= term
    return gf_mul(prod, total)


def _parity_columns_for_data(data_columns: int, redundancy_pct: float) -> int:
    data_columns = max(1, int(data_columns))
    r = max(0.0, min(95.0, float(redundancy_pct))) / 100.0
    if r <= 0:
        return 0
    # Toolkit's config defines redundancy as percentage of each final codeword.
    # Therefore p / (k + p) = r  =>  p = k*r/(1-r).
    return max(1, int(math.ceil(data_columns * r / max(1e-12, 1.0 - r))))


@dataclass
class ToolkitEncodeResult:
    dna: str
    strands: List[str]
    meta: Dict[str, Any]


def encode_toolkit_rs(
    data: bytes,
    *,
    redundancy_pct: float = 15.0,
    rs_rows: int = RS_ROWS,
) -> ToolkitEncodeResult:
    if int(rs_rows) != RS_ROWS:
        raise ValueError("This built-in Toolkit baseline currently follows RS_ROWS=14")
    data = bytes(data or b"")
    symbols = _bytes_to_symbols(data)
    if not symbols:
        symbols = [0]
    data_columns = int(math.ceil(len(symbols) / RS_ROWS))
    parity_columns = _parity_columns_for_data(data_columns, redundancy_pct)
    total_columns = data_columns + parity_columns
    if total_columns >= GF_SIZE:
        raise ValueError("Too many Toolkit columns for 16-bit index")

    # Fill RS matrix by column, exactly like Toolkit mapping_scheme=0.
    rows: List[List[int]] = [[0] * data_columns for _ in range(RS_ROWS)]
    for i, sym in enumerate(symbols):
        row = i % RS_ROWS
        col = i // RS_ROWS
        rows[row][col] = int(sym) & GF_MASK

    all_xs = _column_points(total_columns)
    data_xs = all_xs[:data_columns]
    data_weights = _barycentric_weights(data_xs) if data_columns > 1 else [1]

    encoded_rows: List[List[int]] = []
    for row_vals in rows:
        code_vals = list(row_vals)
        if parity_columns:
            for x in all_xs[data_columns:]:
                if data_columns == 1:
                    code_vals.append(row_vals[0])
                else:
                    code_vals.append(_lagrange_eval(data_xs, row_vals, x, data_weights))
        encoded_rows.append(code_vals)

    strands: List[str] = []
    for col in range(total_columns):
        s = [_index_to_dna(col)]
        for row in range(RS_ROWS):
            s.append(_symbol_to_dna(encoded_rows[row][col]))
        strands.append("".join(s))

    dna = "".join(strands)
    meta = {
        "mapping": "Toolkit RS Baseline",
        "mode": "TOOLKIT_RS_BASELINE",
        "rs_rows": RS_ROWS,
        "symbol_size_bits": 16,
        "file_size": len(data),
        "bytes_len": len(data),
        "input_symbols": len(symbols),
        "data_columns": data_columns,
        "parity_columns": parity_columns,
        "total_columns": total_columns,
        "redundancy_pct_final": float(redundancy_pct),
        "parity_over_data_pct": (100.0 * parity_columns / data_columns) if data_columns else 0.0,
        "strand_length_nt": STRAND_NT,
        "index_nt_per_strand": INDEX_NT,
        "payload_nt_per_strand": RS_ROWS * SYMBOL_NT,
        "dna_length_nt": len(dna),
        "rs_parity_nt": parity_columns * RS_ROWS * SYMBOL_NT,
        "index_nt_total": total_columns * INDEX_NT,
        "extra_full_strand_nt": parity_columns * STRAND_NT,
    }
    return ToolkitEncodeResult(dna=dna, strands=strands, meta=meta)


def _parse_toolkit_strands(dna_text: str) -> List[str]:
    text = str(dna_text or "")
    # Prefer line records if the user provided EncodedStrands/ReconstructedStrands style text.
    records: List[str] = []
    for line in text.splitlines():
        seq = clean_dna(line)
        if seq:
            if len(seq) >= STRAND_NT:
                # Official strands are fixed length. If a line is longer because of accidental
                # whitespace/concatenation, split into chunks.
                for i in range(0, len(seq), STRAND_NT):
                    chunk = seq[i:i + STRAND_NT]
                    if len(chunk) == STRAND_NT:
                        records.append(chunk)
            else:
                # Ignore partial lines.
                pass
    if records:
        return records
    seq = clean_dna(text)
    return [seq[i:i + STRAND_NT] for i in range(0, len(seq), STRAND_NT) if len(seq[i:i + STRAND_NT]) == STRAND_NT]


def strands_to_text(strands: Sequence[str]) -> str:
    return "\n".join(clean_dna(s) for s in strands if clean_dna(s)) + ("\n" if strands else "")


def decode_toolkit_rs(
    dna_text: str,
    *,
    file_size: Optional[int] = None,
    data_columns: Optional[int] = None,
    parity_columns: Optional[int] = None,
    redundancy_pct: float = 15.0,
) -> Tuple[bytes, Dict[str, Any]]:
    strands = _parse_toolkit_strands(dna_text)
    if not strands:
        raise ValueError("No valid 120-nt Toolkit strands were found")

    by_col: Dict[int, List[int]] = {}
    duplicate_columns = 0
    for strand in strands:
        try:
            col = _dna_to_index(strand[:INDEX_NT])
            vals = []
            for row in range(RS_ROWS):
                start = INDEX_NT + row * SYMBOL_NT
                vals.append(_dna_to_symbol(strand[start:start + SYMBOL_NT]))
            if col in by_col:
                duplicate_columns += 1
                # Keep the first clean/reconstructed strand for now.
                continue
            by_col[col] = vals
        except Exception:
            continue

    if not by_col:
        raise ValueError("Toolkit strands could not be parsed into indexed columns")

    max_col = max(by_col)
    if data_columns is None:
        if file_size is not None:
            data_symbols = int(math.ceil(max(0, int(file_size)) / 2.0)) or 1
            data_columns = int(math.ceil(data_symbols / RS_ROWS))
        elif parity_columns is not None:
            data_columns = max(1, max_col + 1 - int(parity_columns))
        else:
            # Last resort: treat all observed columns as data. This cannot remove RS parity
            # and is only useful for skip-ECC/import debugging.
            data_columns = max_col + 1
    data_columns = int(data_columns)
    if parity_columns is None:
        parity_columns = _parity_columns_for_data(data_columns, redundancy_pct)
    parity_columns = int(parity_columns)
    total_columns = data_columns + parity_columns

    # Use only columns inside the expected codeword.
    known_cols = sorted(c for c in by_col if 0 <= c < total_columns)
    if len(known_cols) < data_columns:
        raise ValueError(
            f"Not enough Toolkit columns to decode: have {len(known_cols)}, need {data_columns}"
        )

    all_xs = _column_points(total_columns)
    data_xs = all_xs[:data_columns]
    known_xs = [all_xs[c] for c in known_cols]
    known_weights = _barycentric_weights(known_xs) if len(known_xs) > 1 else [1]

    decoded_rows: List[List[int]] = [[0] * data_columns for _ in range(RS_ROWS)]
    erasure_columns = [c for c in range(total_columns) if c not in by_col]
    missing_data_columns = [c for c in range(data_columns) if c not in by_col]

    for row in range(RS_ROWS):
        known_ys = [by_col[c][row] for c in known_cols]
        for col in range(data_columns):
            if col in by_col:
                decoded_rows[row][col] = by_col[col][row]
            else:
                if len(known_xs) == 1:
                    decoded_rows[row][col] = known_ys[0]
                else:
                    decoded_rows[row][col] = _lagrange_eval(known_xs, known_ys, data_xs[col], known_weights)

    # Convert matrix back by column, matching Toolkit's WriteRSMatrixByColumn.
    data_symbols_count = (int(math.ceil(int(file_size) / 2.0)) if file_size is not None else data_columns * RS_ROWS)
    symbols: List[int] = []
    for i in range(data_symbols_count):
        row = i % RS_ROWS
        col = i // RS_ROWS
        if col < data_columns:
            symbols.append(decoded_rows[row][col])

    out_file_size = int(file_size) if file_size is not None else -1
    data = _symbols_to_bytes(symbols, out_file_size)
    meta: Dict[str, Any] = {
        "mapping": "Toolkit RS Baseline",
        "mode": "TOOLKIT_RS_BASELINE",
        "rs_rows": RS_ROWS,
        "symbol_size_bits": 16,
        "file_size": out_file_size if out_file_size >= 0 else len(data),
        "bytes_len": len(data),
        "data_columns": data_columns,
        "parity_columns": parity_columns,
        "total_columns": total_columns,
        "observed_columns": len(by_col),
        "duplicate_columns": duplicate_columns,
        "erasure_columns": len(erasure_columns),
        "missing_data_columns": len(missing_data_columns),
        "corrected_erasure_columns": len(missing_data_columns),
        "strand_length_nt": STRAND_NT,
    }
    return data, meta
