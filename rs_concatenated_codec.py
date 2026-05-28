from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

from gf_rs_generic import GF2m
from utils_bits_v2 import bytes_to_bitstring
from fragments import clean_dna, gc_fraction, max_homopolymer

# dna_rs_coding-style fields:
# Inner RS over GF(2^6): 1 symbol = 6 bits = 3 nt.
# Outer RS over GF(2^14): across strands.
GF6 = GF2m(6, 0x43)       # x^6 + x + 1
GF14 = GF2m(14, 0x402B)   # primitive polynomial for GF(2^14)

BASES = "ACGT"
BASE2VAL = {b: i for i, b in enumerate(BASES)}

# Full RS strand = FBR + Index + Payload + RS parity + RBR = 125 nt
FBR = "ACACGACGCTCTTCCGATCT"      # 20 nt
RBR = "AGATCGGAAGAGCACACGTCT"     # 21 nt

INNER_INDEX_SYMBOLS = 4           # 12 nt
INNER_PAYLOAD_SYMBOLS = 16        # 48 nt
INNER_INFO_SYMBOLS = INNER_INDEX_SYMBOLS + INNER_PAYLOAD_SYMBOLS
INNER_PARITY_SYMBOLS = 8          # 24 nt
INNER_CODE_SYMBOLS = INNER_INFO_SYMBOLS + INNER_PARITY_SYMBOLS
CORE_DNA_LEN = INNER_CODE_SYMBOLS * 3
FULL_STRAND_LEN = len(FBR) + CORE_DNA_LEN + len(RBR)

OUTER_DATA_SEQS = 10
OUTER_PARITY_SEQS = 2
OUTER_TOTAL_SEQS = OUTER_DATA_SEQS + OUTER_PARITY_SEQS
OUTER_PAYLOAD_SYMBOLS = 6         # 6 * 14 = 84 bits = 14 GF(2^6) symbols, padded to 16

# 24-bit reversible whitening mask for the 4-symbol index.
# This prevents index 0 from appearing as AAAAAAAAAAAA.
INDEX_MASK = 0x1B2D3C


@dataclass
class RSEncodeResult:
    dna: str               # coding-level DNA only: Index + Payload + RS parity
    bits: str
    meta: Dict[str, Any]
    strand_rows: List[Dict[str, Any]]


def _bytes_to_bits(data: bytes) -> List[int]:
    bits: List[int] = []
    for b in bytes(data or b""):
        for shift in range(7, -1, -1):
            bits.append((b >> shift) & 1)
    return bits


def _bits_to_bytes(bits: List[int], nbytes: int | None = None) -> bytes:
    bits = list(bits)
    if len(bits) % 8:
        bits += [0] * (8 - (len(bits) % 8))
    out = bytearray()
    for i in range(0, len(bits), 8):
        v = 0
        for bit in bits[i:i + 8]:
            v = (v << 1) | (bit & 1)
        out.append(v)
    return bytes(out[:nbytes] if nbytes is not None else out)


def _bits_to_symbols(bits: List[int], width: int) -> List[int]:
    bits = list(bits)
    if len(bits) % width:
        bits += [0] * (width - (len(bits) % width))
    out: List[int] = []
    for i in range(0, len(bits), width):
        v = 0
        for bit in bits[i:i + width]:
            v = (v << 1) | (bit & 1)
        out.append(v)
    return out


def _symbols_to_bits(symbols: List[int], width: int) -> List[int]:
    bits: List[int] = []
    for s in symbols:
        for shift in range(width - 1, -1, -1):
            bits.append((int(s) >> shift) & 1)
    return bits


def _symbol6_to_dna(s: int) -> str:
    s = int(s) & 0x3F
    return "".join(BASES[(s >> shift) & 0b11] for shift in (4, 2, 0))


def _dna_to_symbol6(tri: str) -> int:
    tri = clean_dna(tri)
    if len(tri) != 3:
        raise ValueError("GF(2^6) symbol must be exactly 3 nt.")
    v = 0
    for ch in tri:
        v = (v << 2) | BASE2VAL[ch]
    return v & 0x3F


def _symbols6_to_dna(symbols: List[int]) -> str:
    return "".join(_symbol6_to_dna(s) for s in symbols)


def _dna_to_symbols6(dna: str) -> List[int]:
    dna = clean_dna(dna)
    if len(dna) % 3 != 0:
        raise ValueError("DNA length must be a multiple of 3 for GF(2^6) symbols.")
    return [_dna_to_symbol6(dna[i:i + 3]) for i in range(0, len(dna), 3)]


def _index_symbols(global_index: int) -> List[int]:
    masked = (int(global_index) ^ INDEX_MASK) & 0xFFFFFF
    bits = [(masked >> shift) & 1 for shift in range(23, -1, -1)]
    return _bits_to_symbols(bits, 6)[:INNER_INDEX_SYMBOLS]


def _symbols_to_index(symbols: List[int]) -> int:
    bits = _symbols_to_bits(symbols[:INNER_INDEX_SYMBOLS], 6)
    masked = 0
    for bit in bits[:24]:
        masked = (masked << 1) | bit
    return (masked ^ INDEX_MASK) & 0xFFFFFF


def _outer_encode_group(data_seqs: List[List[int]]) -> List[List[int]]:
    xs = list(range(1, OUTER_DATA_SEQS + 1))
    all_seqs = [list(s) for s in data_seqs]
    for p in range(OUTER_PARITY_SEQS):
        x = OUTER_DATA_SEQS + 1 + p
        seq: List[int] = []
        for j in range(OUTER_PAYLOAD_SYMBOLS):
            ys = [data_seqs[i][j] for i in range(OUTER_DATA_SEQS)]
            seq.append(GF14.lagrange_eval(xs, ys, x))
        all_seqs.append(seq)
    return all_seqs


def _outer_recover_group(seq_payloads: Dict[int, List[int]]) -> Tuple[List[List[int]], bool]:
    available = [
        (sid + 1, payload)
        for sid, payload in seq_payloads.items()
        if 0 <= sid < OUTER_TOTAL_SEQS and len(payload) == OUTER_PAYLOAD_SYMBOLS
    ]
    if len(available) < OUTER_DATA_SEQS:
        return [], False

    available = available[:OUTER_DATA_SEQS]
    xs = [x for x, _ in available]
    out: List[List[int]] = []
    for data_id in range(OUTER_DATA_SEQS):
        if data_id in seq_payloads and len(seq_payloads[data_id]) == OUTER_PAYLOAD_SYMBOLS:
            out.append(seq_payloads[data_id])
            continue
        target_x = data_id + 1
        seq: List[int] = []
        for j in range(OUTER_PAYLOAD_SYMBOLS):
            ys = [payload[j] for _, payload in available]
            seq.append(GF14.lagrange_eval(xs, ys, target_x))
        out.append(seq)
    return out, True


def _align_to_template(read: str, template: str) -> str:
    read = clean_dna(read)
    template = clean_dna(template)
    if len(read) == len(template):
        return read

    n, m = len(template), len(read)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    bt = [[""] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = i
        bt[i][0] = "D"
    for j in range(1, m + 1):
        dp[0][j] = j
        bt[0][j] = "I"

    for i in range(1, n + 1):
        ti = template[i - 1]
        for j in range(1, m + 1):
            rj = read[j - 1]
            sub = dp[i - 1][j - 1] + (0 if ti == rj else 1)
            dele = dp[i - 1][j] + 1
            ins = dp[i][j - 1] + 1
            best = min(sub, dele, ins)
            dp[i][j] = best
            bt[i][j] = "M" if best == sub else ("D" if best == dele else "I")

    out = []
    i, j = n, m
    while i > 0 or j > 0:
        b = bt[i][j]
        if b == "M":
            out.append(read[j - 1])
            i -= 1
            j -= 1
        elif b == "D":
            out.append("A")
            i -= 1
        else:
            j -= 1

    aligned = "".join(reversed(out))
    if len(aligned) < len(template):
        aligned += "A" * (len(template) - len(aligned))
    return aligned[:len(template)]


def core_from_full_strand(seq: str) -> str:
    seq = clean_dna(seq)
    fbr = clean_dna(FBR)
    rbr = clean_dna(RBR)
    if len(seq) == FULL_STRAND_LEN and seq.startswith(fbr) and seq.endswith(rbr):
        return seq[len(fbr):len(seq) - len(rbr)]
    if len(seq) == FULL_STRAND_LEN:
        return seq[len(fbr):len(seq) - len(rbr)]
    return seq


def full_from_core(core: str) -> str:
    return clean_dna(FBR) + clean_dna(core) + clean_dna(RBR)


def _make_row(global_no: int, gid: int, sid: int, code_symbols: List[int]) -> Dict[str, Any]:
    core = _symbols6_to_dna(code_symbols)
    index_dna = _symbols6_to_dna(code_symbols[:INNER_INDEX_SYMBOLS])
    payload_dna = _symbols6_to_dna(code_symbols[INNER_INDEX_SYMBOLS:INNER_INFO_SYMBOLS])
    parity_dna = _symbols6_to_dna(code_symbols[INNER_INFO_SYMBOLS:])
    full = full_from_core(core)

    return {
        "No.": str(global_no),
        "Type": "RS data" if sid < OUTER_DATA_SEQS else "RS parity",
        "RS design": "Concatenated Reed-Solomon",
        "RS direct": "true",
        "Toolkit RS direct": "true",
        "Group": str(gid),
        "Shard": str(sid),
        "FBR": clean_dna(FBR),
        "Strand index": index_dna,
        "Index": index_dna,
        "Payload": payload_dna,
        "Filler": "",
        "RS parity": parity_dna,
        "RBR": clean_dna(RBR),
        "Core strand": core,
        "Full strand": full,
        "Wet-lab reference payload": core,
        "FBR length": str(len(clean_dna(FBR))),
        "Index length": str(len(index_dna)),
        "Payload length": str(len(payload_dna)),
        "RS parity length": str(len(parity_dna)),
        "Filler length": "0",
        "RBR length": str(len(clean_dna(RBR))),
        "Target total length": str(FULL_STRAND_LEN),
        "Total length": str(len(full)),
        "GC content": f"{gc_fraction(full):.3f}",
        "Longest homopolymer": str(max_homopolymer(full)),
        "DNA": full,
    }


def encode(data: bytes) -> RSEncodeResult:
    raw = bytes(data or b"")
    bits = bytes_to_bitstring(raw)
    blob = len(raw).to_bytes(4, "big") + raw
    outer_symbols = _bits_to_symbols(_bytes_to_bits(blob), 14)

    symbols_per_group = OUTER_DATA_SEQS * OUTER_PAYLOAD_SYMBOLS
    groups = max(1, math.ceil(len(outer_symbols) / symbols_per_group))
    outer_symbols += [0] * (groups * symbols_per_group - len(outer_symbols))

    rows: List[Dict[str, Any]] = []
    core_parts: List[str] = []

    for gid in range(groups):
        gsym = outer_symbols[gid * symbols_per_group:(gid + 1) * symbols_per_group]
        data_seqs = [
            gsym[i * OUTER_PAYLOAD_SYMBOLS:(i + 1) * OUTER_PAYLOAD_SYMBOLS]
            for i in range(OUTER_DATA_SEQS)
        ]
        all_seqs = _outer_encode_group(data_seqs)

        for sid, outer_payload in enumerate(all_seqs):
            payload_bits = _symbols_to_bits(outer_payload, 14)
            payload_inner = _bits_to_symbols(payload_bits, 6)

            # 6 GF(2^14) symbols = 84 bits = 14 GF(2^6) symbols.
            # Pad to 16 symbols so Payload is exactly 48 nt.
            if len(payload_inner) < INNER_PAYLOAD_SYMBOLS:
                payload_inner += [0] * (INNER_PAYLOAD_SYMBOLS - len(payload_inner))
            payload_inner = payload_inner[:INNER_PAYLOAD_SYMBOLS]

            global_index = gid * OUTER_TOTAL_SEQS + sid
            info_symbols = _index_symbols(global_index) + payload_inner
            code_symbols = GF6.rs_encode_msg(info_symbols, INNER_PARITY_SYMBOLS)
            row = _make_row(len(rows) + 1, gid, sid, code_symbols)
            rows.append(row)
            core_parts.append(row["Core strand"])

    coding_dna = "".join(core_parts)
    meta = {
        "mapping": "Reed-Solomon",
        "mode": "REED_SOLOMON",
        "method": "Concatenated Reed-Solomon",
        "bytes_len": len(raw),
        "bits_len": len(bits),
        "data_length": len(raw),
        "groups": groups,
        "inner_field": "GF(2^6)",
        "outer_field": "GF(2^14)",
        "inner_index_symbols": INNER_INDEX_SYMBOLS,
        "inner_payload_symbols": INNER_PAYLOAD_SYMBOLS,
        "inner_info_symbols": INNER_INFO_SYMBOLS,
        "inner_parity_symbols": INNER_PARITY_SYMBOLS,
        "inner_code_symbols": INNER_CODE_SYMBOLS,
        "outer_data_seqs": OUTER_DATA_SEQS,
        "outer_parity_seqs": OUTER_PARITY_SEQS,
        "outer_total_seqs": OUTER_TOTAL_SEQS,
        "outer_payload_symbols": OUTER_PAYLOAD_SYMBOLS,
        "index_mask": INDEX_MASK,
        "fbr": clean_dna(FBR),
        "rbr": clean_dna(RBR),
        "core_strand_length": CORE_DNA_LEN,
        "strand_length": FULL_STRAND_LEN,
        "strand_count": len(rows),
        "rs_strands": rows,
    }
    return RSEncodeResult(dna=coding_dna, bits=bits, meta=meta, strand_rows=rows)


def _decode_one_core(core_dna: str, template_core: str | None = None) -> Tuple[int, int, List[int]]:
    core_dna = clean_dna(core_dna)
    if template_core and len(core_dna) != len(clean_dna(template_core)):
        core_dna = _align_to_template(core_dna, template_core)

    symbols = _dna_to_symbols6(core_dna)
    info_symbols = GF6.rs_correct_msg(symbols, INNER_PARITY_SYMBOLS)
    global_index = _symbols_to_index(info_symbols[:INNER_INDEX_SYMBOLS])
    gid = global_index // OUTER_TOTAL_SEQS
    sid = global_index % OUTER_TOTAL_SEQS
    payload_inner = info_symbols[INNER_INDEX_SYMBOLS:INNER_INFO_SYMBOLS]

    # Only first 84 bits encode the 6 outer GF(2^14) symbols.
    payload_bits = _symbols_to_bits(payload_inner, 6)
    outer_payload = _bits_to_symbols(payload_bits[:OUTER_PAYLOAD_SYMBOLS * 14], 14)
    return gid, sid, outer_payload


def decode(dna: str, codec_meta: Dict[str, Any] | None = None) -> Tuple[bytes, Dict[str, Any]]:
    codec_meta = codec_meta or {}
    data_len = int(codec_meta.get("data_length", codec_meta.get("bytes_len", 0)))
    groups = int(codec_meta.get("groups", 0))

    dna = clean_dna(dna)
    core_len = int(codec_meta.get("core_strand_length", CORE_DNA_LEN))
    full_len = int(codec_meta.get("strand_length", FULL_STRAND_LEN))

    chunks: List[str] = []
    if len(dna) >= full_len and len(dna) % full_len == 0:
        for i in range(0, len(dna), full_len):
            chunks.append(core_from_full_strand(dna[i:i + full_len]))
    else:
        for i in range(0, len(dna), core_len):
            c = dna[i:i + core_len]
            if len(c) == core_len:
                chunks.append(c)

    if groups <= 0:
        groups = max(1, math.ceil(len(chunks) / OUTER_TOTAL_SEQS))

    seq_payloads_by_group: Dict[int, Dict[int, List[int]]] = {gid: {} for gid in range(groups)}
    report_rows: List[Dict[str, Any]] = []
    corrected = 0
    failed = 0

    for i, chunk in enumerate(chunks, start=1):
        status = "Fail"
        gid: Any = ""
        sid: Any = ""
        try:
            gid, sid, outer_payload = _decode_one_core(chunk)
            if 0 <= gid < groups and 0 <= sid < OUTER_TOTAL_SEQS:
                seq_payloads_by_group[gid][sid] = outer_payload
                status = "Success"
                corrected += 1
            else:
                failed += 1
        except Exception:
            failed += 1
        report_rows.append({"No.": i, "Group": gid, "Shard": sid, "Status": status})

    recovered_outer_symbols: List[int] = []
    recovered_groups = 0
    failed_groups = 0

    for gid in range(groups):
        data_seqs, ok = _outer_recover_group(seq_payloads_by_group.get(gid, {}))
        if not ok:
            failed_groups += 1
            continue
        recovered_groups += 1
        for seq in data_seqs:
            recovered_outer_symbols.extend(seq)

    bits = _symbols_to_bits(recovered_outer_symbols, 14)
    blob = _bits_to_bytes(bits) if bits else b""

    if len(blob) >= 4:
        declared_len = int.from_bytes(blob[:4], "big")
        if declared_len < 0 or declared_len > max(0, len(blob) - 4):
            declared_len = data_len
        data = blob[4:4 + declared_len]
    else:
        data = b""

    meta = {
        "mode": "REED_SOLOMON",
        "method": "Concatenated Reed-Solomon",
        "bytes_len": len(data),
        "bits_len": len(bytes_to_bitstring(data)),
        "strand_chunks_seen": len(chunks),
        "corrected_strands": corrected,
        "failed_strands": failed,
        "recovered_groups": recovered_groups,
        "failed_groups": failed_groups,
        "recovery_table": report_rows,
        "valid": (failed_groups == 0 and (data_len == 0 or len(data) == data_len)),
    }
    return data, meta


def strand_rows_from_meta(codec_meta: Dict[str, Any]) -> List[Dict[str, Any]]:
    return list((codec_meta or {}).get("rs_strands", []) or [])
