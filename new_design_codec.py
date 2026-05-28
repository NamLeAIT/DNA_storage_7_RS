from __future__ import annotations

import math
import zlib
from collections import Counter, defaultdict
from typing import Any, Dict, List, Optional, Tuple

BASES = "ATGC"
VAL2BASE = ["A", "T", "G", "C"]
BASE2VAL = {b: i for i, b in enumerate(VAL2BASE)}

# Same FSM tables as the prototype.
F = [
    [1, 2, 3, 0],
    [0, 1, 3, 2],
    [2, 3, 1, 0],
    [3, 0, 2, 1],
]
G = [[(F[s][x] + s) % 4 for x in range(4)] for s in range(4)]

DEFAULT_DATA_BLOCK_LEN = 70
CRC8_SYMBOLS = 4
STATE_SYMBOLS = 3
FULL_REPAIR_BLOCK_LEN = DEFAULT_DATA_BLOCK_LEN + CRC8_SYMBOLS + STATE_SYMBOLS


def clean_dna(seq: str) -> str:
    return "".join(ch for ch in str(seq).upper() if ch in BASE2VAL)


def vals_to_dna(vals: List[int]) -> str:
    return "".join(VAL2BASE[int(v) & 3] for v in vals)


def dna_to_vals(dna: str) -> List[int]:
    return [BASE2VAL[ch] for ch in clean_dna(dna)]


def flatten_blocks(blocks: List[List[int]]) -> List[int]:
    out: List[int] = []
    for b in blocks:
        out.extend(b)
    return out


def bytes_to_vals(data: bytes) -> List[int]:
    vals: List[int] = []
    for byte in data:
        vals.append((byte >> 6) & 3)
        vals.append((byte >> 4) & 3)
        vals.append((byte >> 2) & 3)
        vals.append(byte & 3)
    return vals


def vals_to_bytes(vals: List[int]) -> bytes:
    vals = [int(v) & 3 for v in vals]
    usable = (len(vals) // 4) * 4
    vals = vals[:usable]
    out = bytearray()
    for i in range(0, len(vals), 4):
        b = (vals[i] << 6) | (vals[i + 1] << 4) | (vals[i + 2] << 2) | vals[i + 3]
        out.append(b)
    return bytes(out)


def crc8_bytes(data: bytes, poly: int = 0x07, init: int = 0x00) -> int:
    crc = init
    for b in data:
        crc ^= b
        for _ in range(8):
            if crc & 0x80:
                crc = ((crc << 1) ^ poly) & 0xFF
            else:
                crc = (crc << 1) & 0xFF
    return crc & 0xFF


def vals_crc8(vals: List[int]) -> int:
    return crc8_bytes(bytes(int(v) & 3 for v in vals))


def int_to_vals_2bit(x: int, bit_len: int) -> List[int]:
    bits = format(int(x), f"0{bit_len}b")
    return [int(bits[i:i + 2], 2) for i in range(0, bit_len, 2)]


def vals_2bit_to_int(vals: List[int]) -> int:
    bits = "".join(f"{int(v) & 3:02b}" for v in vals)
    return int(bits, 2) if bits else 0


def append_block_crc8(block_vals: List[int]) -> Tuple[List[int], int, List[int]]:
    crc = vals_crc8(block_vals)
    crc_vals = int_to_vals_2bit(crc, 8)
    return block_vals + crc_vals, crc, crc_vals


def split_block_payload_crc8(block_vals_with_crc8: List[int]) -> Tuple[List[int], int, List[int]]:
    if len(block_vals_with_crc8) < CRC8_SYMBOLS:
        return [], 0, []
    block_vals = block_vals_with_crc8[:-CRC8_SYMBOLS]
    crc_vals = block_vals_with_crc8[-CRC8_SYMBOLS:]
    stored_crc = vals_2bit_to_int(crc_vals)
    return block_vals, stored_crc, crc_vals


def embedded_block_crc8_ok(block_vals_with_crc8: List[int]) -> bool:
    block_vals, stored_crc, _ = split_block_payload_crc8(block_vals_with_crc8)
    return vals_crc8(block_vals) == stored_crc


def encode_once(payload: List[int], s0: int = 0) -> Tuple[List[int], int]:
    s = int(s0) & 3
    out: List[int] = []
    for x in payload:
        x = int(x) & 3
        y = F[s][x]
        out.append(y)
        s = G[s][x]
    return out, s


def encode_3layer(payload: List[int]) -> Tuple[List[int], List[int], List[int], Tuple[int, int, int]]:
    dna1, s1 = encode_once(payload)
    dna2, s2 = encode_once(dna1)
    dna3, s3 = encode_once(dna2)
    return dna1, dna2, dna3, (s1, s2, s3)


def composite_step(state: Tuple[int, int, int], x: int) -> Tuple[Tuple[int, int, int], int]:
    s1, s2, s3 = state
    x = int(x) & 3
    y1 = F[s1][x]
    ns1 = G[s1][x]
    y2 = F[s2][y1]
    ns2 = G[s2][y1]
    y3 = F[s3][y2]
    ns3 = G[s3][y2]
    return (ns1, ns2, ns3), y3


def suffix_feasibility(dna_obs: List[int], target_state: Tuple[int, int, int], max_errors: int):
    n = len(dna_obs)
    all_states = [(a, b, c) for a in range(4) for b in range(4) for c in range(4)]
    feas = [defaultdict(Counter) for _ in range(n + 1)]
    feas[n][target_state][0] = 1
    for i in range(n - 1, -1, -1):
        obs = dna_obs[i]
        for state in all_states:
            ctr = Counter()
            for x in range(4):
                next_state, y3 = composite_step(state, x)
                step_cost = 0 if y3 == obs else 1
                for rem_cost, count in feas[i + 1].get(next_state, {}).items():
                    total = step_cost + rem_cost
                    if total <= max_errors:
                        ctr[total] += count
            if ctr:
                feas[i][state] = ctr
    return feas


def enumerate_block_candidates(
    dna_obs: List[int],
    target_state: Tuple[int, int, int],
    max_errors: int,
    require_crc8: bool = True,
    max_candidates: int = 64,
) -> List[Dict[str, Any]]:
    n = len(dna_obs)
    start_state = (0, 0, 0)
    feas = suffix_feasibility(dna_obs, target_state, max_errors)
    candidates: List[Dict[str, Any]] = []

    def dfs(i: int, state: Tuple[int, int, int], cost: int, block_prefix: List[int], dna3_prefix: List[int]):
        if cost > max_errors or len(candidates) >= max_candidates:
            return
        if i == n:
            if state != target_state:
                return
            if require_crc8 and not embedded_block_crc8_ok(block_prefix):
                return
            block_data, stored_crc8, crc8_vals = split_block_payload_crc8(block_prefix)
            candidates.append({
                "block_data": block_data,
                "block_with_crc8": block_prefix.copy(),
                "dna3": dna3_prefix.copy(),
                "distance": cost,
                "stored_crc8": stored_crc8,
                "crc8_vals": crc8_vals,
            })
            return

        obs = dna_obs[i]
        for x in range(4):
            next_state, y3 = composite_step(state, x)
            new_cost = cost + (0 if y3 == obs else 1)
            if new_cost > max_errors:
                continue
            possible = False
            for rem_cost in feas[i + 1].get(next_state, {}):
                if new_cost + rem_cost <= max_errors:
                    possible = True
                    break
            if not possible:
                continue
            block_prefix.append(x)
            dna3_prefix.append(y3)
            dfs(i + 1, next_state, new_cost, block_prefix, dna3_prefix)
            block_prefix.pop()
            dna3_prefix.pop()

    dfs(0, start_state, 0, [], [])
    candidates.sort(key=lambda c: (c["distance"], c["block_data"]))
    return candidates


def _states_by_distance(obs_state: Tuple[int, int, int], max_dist: int) -> List[Tuple[Tuple[int, int, int], int]]:
    states = [(a, b, c) for a in range(4) for b in range(4) for c in range(4)]
    out = []
    for st in states:
        d = sum(1 for a, b in zip(obs_state, st) if a != b)
        if d <= max_dist:
            out.append((st, d))
    out.sort(key=lambda item: (item[1], item[0]))
    return out


def encode_new_design(data: bytes, block_data_len: int = DEFAULT_DATA_BLOCK_LEN) -> Dict[str, Any]:
    vals = bytes_to_vals(data)
    blocks = [vals[i:i + block_data_len] for i in range(0, len(vals), block_data_len)]
    repair_blocks: List[str] = []
    block_reports: List[Dict[str, Any]] = []

    for i, block in enumerate(blocks, start=1):
        block_with_crc8, crc8_value, crc8_vals = append_block_crc8(block)
        _, _, dna3_vals, target_state = encode_3layer(block_with_crc8)
        state_vals = list(target_state)
        repair_vals = dna3_vals + state_vals
        repair_blocks.append(vals_to_dna(repair_vals))
        block_reports.append({
            "block": i,
            "data_symbols": len(block),
            "crc8_symbols": CRC8_SYMBOLS,
            "state_symbols": STATE_SYMBOLS,
            "repair_symbols": len(repair_vals),
            "crc8": crc8_value,
            "state": state_vals,
        })

    dna = "".join(repair_blocks)
    overhead = len(dna) - len(vals)
    meta = {
        "design": "New Design",
        "block_data_len": int(block_data_len),
        "crc8_symbols": CRC8_SYMBOLS,
        "state_symbols": STATE_SYMBOLS,
        "full_repair_block_len": int(block_data_len) + CRC8_SYMBOLS + STATE_SYMBOLS,
        "original_bytes": len(data),
        "payload_symbols": len(vals),
        "num_blocks": len(blocks),
        "repair_symbols": len(dna),
        "overhead_symbols": overhead,
        "expansion_factor": (len(dna) / len(vals)) if vals else 1.0,
        "state_embedded": True,
    }
    return {"dna": dna, "meta": meta, "blocks": block_reports}


def _split_repair_blocks(dna: str, block_data_len: int = DEFAULT_DATA_BLOCK_LEN) -> List[List[int]]:
    vals = dna_to_vals(dna)
    full_len = int(block_data_len) + CRC8_SYMBOLS + STATE_SYMBOLS
    blocks: List[List[int]] = []
    pos = 0
    while pos < len(vals):
        remaining = len(vals) - pos
        take = min(full_len, remaining)
        block = vals[pos:pos + take]
        if len(block) < (CRC8_SYMBOLS + STATE_SYMBOLS + 1):
            # Ignore incomplete trailing filler-like sequence.
            break
        blocks.append(block)
        pos += take
    return blocks


def _decode_one_repair_block(
    repair_block_vals: List[int],
    per_block_max_errors: int,
    max_candidates: int = 128,
) -> List[Dict[str, Any]]:
    if len(repair_block_vals) < CRC8_SYMBOLS + STATE_SYMBOLS + 1:
        return []
    dna_obs = repair_block_vals[:-STATE_SYMBOLS]
    obs_state_vals = repair_block_vals[-STATE_SYMBOLS:]
    obs_state = tuple(int(v) & 3 for v in obs_state_vals)  # type: ignore[assignment]
    out: List[Dict[str, Any]] = []

    for target_state, state_dist in _states_by_distance(obs_state, per_block_max_errors):
        remaining = per_block_max_errors - state_dist
        if remaining < 0:
            continue
        cands = enumerate_block_candidates(
            dna_obs=dna_obs,
            target_state=target_state,
            max_errors=remaining,
            require_crc8=True,
            max_candidates=max_candidates,
        )
        for c in cands:
            total_dist = int(c["distance"]) + state_dist
            out.append({
                **c,
                "distance": total_dist,
                "fsm_distance": int(c["distance"]),
                "state_distance": state_dist,
                "target_state": list(target_state),
                "observed_state": list(obs_state),
            })
        if out and state_dist == 0:
            # Prefer the embedded state when it works; this keeps decoding fast.
            break

    out.sort(key=lambda c: (c["distance"], c.get("state_distance", 0), c["block_data"]))
    return out[:max_candidates]


def decode_new_design(
    dna: str,
    block_data_len: int = DEFAULT_DATA_BLOCK_LEN,
    per_block_schedule: Tuple[int, ...] = (0, 1, 2, 3),
    candidate_limit: int = 24,
) -> Dict[str, Any]:
    repair_blocks = _split_repair_blocks(dna, block_data_len=block_data_len)
    if not repair_blocks:
        raise ValueError("No valid New Design repair blocks found.")

    block_candidate_sets: List[List[Dict[str, Any]]] = []
    used_per_block = None
    for max_err in per_block_schedule:
        sets = []
        ok = True
        for block_vals in repair_blocks:
            cands = _decode_one_repair_block(block_vals, per_block_max_errors=max_err, max_candidates=128)
            if not cands:
                ok = False
                break
            sets.append(cands)
        if ok:
            block_candidate_sets = sets
            used_per_block = max_err
            break

    if not block_candidate_sets:
        raise ValueError("New Design repair failed: at least one block could not pass block check.")

    # Beam-combine block candidates into a small possible-output library.
    beam: List[Dict[str, Any]] = [{"vals": [], "distance": 0, "block_choices": []}]
    for block_idx, cands in enumerate(block_candidate_sets, start=1):
        new_beam: List[Dict[str, Any]] = []
        for prefix in beam:
            for cand in cands[: min(len(cands), max(candidate_limit, 24))]:
                new_beam.append({
                    "vals": prefix["vals"] + cand["block_data"],
                    "distance": prefix["distance"] + int(cand["distance"]),
                    "block_choices": prefix["block_choices"] + [{
                        "block": block_idx,
                        "distance": int(cand["distance"]),
                        "state_distance": int(cand.get("state_distance", 0)),
                    }],
                })
        new_beam.sort(key=lambda c: (c["distance"], len(c["vals"])))
        beam = new_beam[:candidate_limit]

    candidates: List[Dict[str, Any]] = []
    for rank, cand in enumerate(beam, start=1):
        data = vals_to_bytes(cand["vals"])
        candidates.append({
            "rank": rank,
            "data": data,
            "distance": int(cand["distance"]),
            "bytes_len": len(data),
            "vals_len": len(cand["vals"]),
            "block_choices": cand["block_choices"],
        })

    best = candidates[0]
    return {
        "status": "best_output" if best["distance"] else "recovered",
        "best_data": best["data"],
        "best_distance": best["distance"],
        "corrected_bases": best["distance"],
        "num_blocks": len(repair_blocks),
        "used_per_block_max": used_per_block,
        "candidate_count": len(candidates),
        "candidates": candidates,
        "meta": {
            "design": "New Design",
            "block_data_len": int(block_data_len),
            "num_blocks": len(repair_blocks),
            "used_per_block_max": used_per_block,
            "state_embedded": True,
        },
    }
