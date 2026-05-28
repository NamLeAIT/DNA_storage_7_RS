from __future__ import annotations

import csv
import hashlib
import io
import math
import random
from typing import Dict, List, Set

import dna_codec


BASES = "ACGT"


def clean_dna(seq: str) -> str:
    return dna_codec.clean_dna_text(seq)


def split_dna(seq: str, payload_len: int) -> List[str]:
    seq = clean_dna(seq)
    if payload_len <= 0:
        payload_len = 100
    return [seq[i:i + payload_len] for i in range(0, len(seq), payload_len)]


def max_homopolymer(seq: str) -> int:
    seq = clean_dna(seq)
    if not seq:
        return 0

    cur = 1
    mx = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur += 1
            mx = max(mx, cur)
        else:
            cur = 1
    return mx


def gc_fraction(seq: str) -> float:
    seq = clean_dna(seq)
    if not seq:
        return 0.0
    return sum(b in "GC" for b in seq) / len(seq)


def _hash_to_dna(seed: str, length: int) -> str:
    """
    Deterministically generate a DNA candidate from a seed.
    This is not random at runtime; same seed gives same sequence.
    """
    digest = hashlib.sha256(seed.encode("utf-8")).digest()
    bits = []
    for byte in digest:
        for shift in (6, 4, 2, 0):
            bits.append((byte >> shift) & 0b11)

    out = []
    i = 0
    while len(out) < length:
        if i >= len(bits):
            digest = hashlib.sha256((seed + "|" + str(i)).encode("utf-8")).digest()
            bits = []
            for byte in digest:
                for shift in (6, 4, 2, 0):
                    bits.append((byte >> shift) & 0b11)
            i = 0

        out.append(BASES[bits[i]])
        i += 1

    return "".join(out)


def _index_score(index_seq: str, left_context: str, right_context: str) -> float:
    """
    Lower score is better.
    Penalizes homopolymers, poor GC balance, and boundary repeats.
    """
    full = clean_dna(left_context) + clean_dna(index_seq) + clean_dna(right_context)
    idx = clean_dna(index_seq)

    score = 0.0

    hp_full = max_homopolymer(full)
    hp_idx = max_homopolymer(idx)

    score += max(0, hp_full - 2) * 100.0
    score += max(0, hp_idx - 2) * 80.0

    gc = gc_fraction(idx)
    score += abs(gc - 0.50) * 20.0

    # Boundary penalty: avoid FBR ending with same base as index start
    if left_context and idx and clean_dna(left_context)[-1:] == idx[:1]:
        score += 10.0

    # Boundary penalty: avoid index ending with same base as payload/RBR start
    if right_context and idx and idx[-1:] == clean_dna(right_context)[:1]:
        score += 10.0

    return score


def generate_strand_index(
    strand_no: int,
    index_len: int,
    used_indices: Set[str],
    left_context: str = "",
    right_context: str = "",
    max_hp_allowed: int = 2,
) -> str:
    """
    Generate a unique DNA strand index with low homopolymer tendency.

    It tries multiple deterministic candidates and selects the best one.
    The index is checked together with neighboring sequence context:
        FBR + strand_index + payload_start
    """
    if index_len <= 0:
        return ""

    best_seq = ""
    best_score = float("inf")

    for salt in range(5000):
        seed = f"strand_index|{strand_no}|{index_len}|{salt}"
        cand = _hash_to_dna(seed, index_len)

        if cand in used_indices:
            continue

        score = _index_score(cand, left_context, right_context)

        if score < best_score:
            best_seq = cand
            best_score = score

        full_context = clean_dna(left_context) + cand + clean_dna(right_context)
        if (
            score <= 1e-9
            and max_homopolymer(cand) <= max_hp_allowed
            and max_homopolymer(full_context) <= max_hp_allowed
        ):
            used_indices.add(cand)
            return cand

    if not best_seq:
        raise RuntimeError("Could not generate a unique strand index.")

    used_indices.add(best_seq)
    return best_seq


def choose_auto_strand_design(
    dna_len: int,
    fbr_len: int,
    rbr_len: int,
    index_len: int,
    min_total_len: int = 120,
    max_total_len: int = 130,
) -> Dict[str, int]:
    """
    Choose a strand total length in [min_total_len, max_total_len].

    Objective:
      1) valid payload capacity > 0
      2) minimize filler nt added at the end
      3) prefer larger payload capacity if filler tie
      4) prefer total length close to the middle of the requested range
    """
    dna_len = max(0, int(dna_len))
    min_total_len = int(min_total_len)
    max_total_len = int(max_total_len)
    if min_total_len > max_total_len:
        min_total_len, max_total_len = max_total_len, min_total_len

    fixed_len = int(fbr_len) + int(index_len) + int(rbr_len)
    mid = (min_total_len + max_total_len) / 2.0
    best = None

    for total_len in range(min_total_len, max_total_len + 1):
        payload_capacity = total_len - fixed_len
        if payload_capacity <= 0:
            continue
        n_strands = max(1, math.ceil(dna_len / payload_capacity)) if dna_len else 0
        filler_total = (n_strands * payload_capacity - dna_len) if n_strands else 0
        score = (filler_total, -payload_capacity, abs(total_len - mid), n_strands)
        if best is None or score < best[0]:
            best = (
                score,
                {
                    "target_total_len": total_len,
                    "payload_capacity": payload_capacity,
                    "n_strands": n_strands,
                    "filler_total": filler_total,
                    "fixed_len": fixed_len,
                },
            )

    if best is None:
        raise ValueError(
            "No valid total strand length. Increase total length or reduce FBR/RBR/index length."
        )
    return best[1]


def _sequence_score(seq: str, left_context: str = "", right_context: str = "") -> float:
    """Lower score is better for filler/index-like helper sequences."""
    seq = clean_dna(seq)
    full = clean_dna(left_context) + seq + clean_dna(right_context)
    if not seq:
        return 0.0

    score = 0.0
    hp_full = max_homopolymer(full)
    hp_seq = max_homopolymer(seq)
    score += max(0, hp_full - 2) * 100.0
    score += max(0, hp_seq - 2) * 80.0
    score += abs(gc_fraction(seq) - 0.50) * 20.0

    if left_context and seq and clean_dna(left_context)[-1:] == seq[:1]:
        score += 10.0
    if right_context and seq and seq[-1:] == clean_dna(right_context)[:1]:
        score += 10.0
    return score


def generate_filler(
    strand_no: int,
    filler_len: int,
    left_context: str = "",
    right_context: str = "",
    max_hp_allowed: int = 2,
) -> str:
    """
    Deterministically generate filler nt.

    Filler is not payload. It is only used to make full strands equal length.
    The sequence is selected to reduce homopolymers and keep GC close to 50%.
    """
    filler_len = int(filler_len)
    if filler_len <= 0:
        return ""

    best_seq = ""
    best_score = float("inf")
    for salt in range(5000):
        cand = _hash_to_dna(f"filler|{strand_no}|{filler_len}|{salt}", filler_len)
        score = _sequence_score(cand, left_context, right_context)
        if score < best_score:
            best_seq = cand
            best_score = score

        full_context = clean_dna(left_context) + cand + clean_dna(right_context)
        if score <= 1e-9 and max_homopolymer(full_context) <= max_hp_allowed:
            return cand

    return best_seq


def prepare_dna_strands(
    dna: str,
    payload_len: int | None = None,
    fbr: str = "",
    rbr: str = "",
    index_len: int = 8,
    target_total_len: int | None = None,
    add_filler: bool = True,
) -> List[Dict[str, str]]:
    """
    Prepare DNA strands as:
        FBR + strand_index + payload + filler + RBR

    Notes:
      - Payload is the real encoded DNA data.
      - Filler is only padding to equalize strand length and must be ignored during decoding.
      - If target_total_len is provided, payload capacity is computed from the total length.
      - If target_total_len is not provided, payload_len is used for backward compatibility.
    """
    dna = clean_dna(dna)
    fbr = clean_dna(fbr)
    rbr = clean_dna(rbr)
    index_len = int(index_len)

    if target_total_len is not None:
        target_total_len = int(target_total_len)
        payload_capacity = target_total_len - len(fbr) - index_len - len(rbr)
        if payload_capacity <= 0:
            raise ValueError(
                "Total strand length is too short for FBR + index + RBR. "
                f"total={target_total_len}, FBR={len(fbr)}, index={index_len}, RBR={len(rbr)}"
            )
    else:
        payload_capacity = int(payload_len or 100)
        target_total_len = len(fbr) + index_len + payload_capacity + len(rbr)

    payloads = split_dna(dna, payload_capacity)
    used_indices: Set[str] = set()
    rows: List[Dict[str, str]] = []

    for i, payload in enumerate(payloads, start=1):
        right_context_for_index = payload[:8] if payload else rbr[:8]

        strand_index = generate_strand_index(
            strand_no=i,
            index_len=index_len,
            used_indices=used_indices,
            left_context=fbr[-8:],
            right_context=right_context_for_index,
            max_hp_allowed=2,
        )

        filler_len = max(0, payload_capacity - len(payload)) if add_filler else 0
        filler = generate_filler(
            strand_no=i,
            filler_len=filler_len,
            left_context=payload[-8:],
            right_context=rbr[:8],
            max_hp_allowed=2,
        )

        full_strand = fbr + strand_index + payload + filler + rbr
        hp_full = dna_codec.homopolymer_stats(full_strand)

        rows.append(
            {
                "No.": str(i),
                "FBR": fbr,
                "Strand index": strand_index,
                "Payload": payload,
                "Filler": filler,
                "RBR": rbr,
                "Full strand": full_strand,
                "FBR length": str(len(fbr)),
                "Index length": str(len(strand_index)),
                "Payload length": str(len(payload)),
                "Filler length": str(len(filler)),
                "RBR length": str(len(rbr)),
                "Payload capacity": str(payload_capacity),
                "Target total length": str(target_total_len),
                "Total length": str(len(full_strand)),
                "Original DNA length": str(len(dna)),
                "GC content": f"{gc_fraction(full_strand):.4f}",
                "Longest homopolymer": str(hp_full.get("longest", 0)),
                "Homopolymer count": str(hp_full.get("count_ge2", 0)),
            }
        )

    return rows


def mutate_dna_sequence(
    seq: str,
    substitution_rate: float = 0.0,
    insertion_rate: float = 0.0,
    deletion_rate: float = 0.0,
    seed: int | str = 1,
    allow_indels: bool = False,
) -> tuple[str, Dict[str, int]]:
    """
    Substitute bases in one DNA sequence.

    Default behavior is substitution only, so sequence length stays unchanged.
    Insertion/deletion are ignored unless allow_indels=True is explicitly passed.
    """
    seq = clean_dna(seq)
    rng = random.Random(str(seed))
    out: List[str] = []
    counts: Dict[str, int] = {
        "Substitute count": 0,
        "Changed bases": 0,
        "Error count": 0,
    }

    for base in seq:
        if allow_indels and rng.random() < float(deletion_rate):
            counts["Deletion count"] = counts.get("Deletion count", 0) + 1
            counts["Error count"] += 1
            if rng.random() < float(insertion_rate):
                out.append(rng.choice(BASES))
                counts["Insertion count"] = counts.get("Insertion count", 0) + 1
                counts["Error count"] += 1
            continue

        if rng.random() < float(substitution_rate):
            choices = [b for b in BASES if b != base]
            out.append(rng.choice(choices))
            counts["Substitute count"] += 1
            counts["Changed bases"] += 1
            counts["Error count"] += 1
        else:
            out.append(base)

        if allow_indels and rng.random() < float(insertion_rate):
            out.append(rng.choice(BASES))
            counts["Insertion count"] = counts.get("Insertion count", 0) + 1
            counts["Error count"] += 1

    return "".join(out), counts


def add_errors_to_strand_rows(
    rows: List[Dict[str, str]],
    substitution_rate: float = 0.0,
    insertion_rate: float = 0.0,
    deletion_rate: float = 0.0,
    seed: int = 1,
    scope: str = "Payload only",
    allow_indels: bool = False,
) -> List[Dict[str, str]]:
    """
    Add substitute changes after strand preparation.

    Default behavior:
      - substitute payload bases only
      - keep sequence length unchanged
      - do not run insertion/deletion unless allow_indels=True is explicitly passed
    """
    out_rows: List[Dict[str, str]] = []
    for row in rows:
        new = dict(row)
        row_no = int(str(row.get("No.", "0") or "0"))
        row_seed = int(seed) + row_no * 1000003
        new["Error scope"] = str(scope)

        if scope == "Full strand":
            mutated, counts = mutate_dna_sequence(
                row.get("Full strand", ""),
                substitution_rate=substitution_rate,
                insertion_rate=insertion_rate,
                deletion_rate=deletion_rate,
                seed=row_seed,
                allow_indels=bool(allow_indels),
            )
            new["Substitute full strand"] = mutated
            new["Substitute payload"] = ""
        else:
            mutated_payload, counts = mutate_dna_sequence(
                row.get("Payload", ""),
                substitution_rate=substitution_rate,
                insertion_rate=insertion_rate,
                deletion_rate=deletion_rate,
                seed=row_seed,
                allow_indels=bool(allow_indels),
            )
            new["Substitute payload"] = mutated_payload
            new["Substitute full strand"] = (
                row.get("FBR", "")
                + row.get("Strand index", "")
                + mutated_payload
                + row.get("Filler", "")
                + row.get("RBR", "")
            )

        # Preserve the old Error count field for compatibility, but show only Substitute in the UI.
        for k, v in counts.items():
            new[k] = str(v)

        if new.get("Substitute payload"):
            ep = new.get("Substitute payload", "")
            new["Substitute payload length"] = str(len(clean_dna(ep)))
        if new.get("Substitute full strand"):
            ef = new.get("Substitute full strand", "")
            hp_ef = dna_codec.homopolymer_stats(ef)
            new["Substitute full length"] = str(len(clean_dna(ef)))
            new["Substitute GC content"] = f"{gc_fraction(ef):.4f}"
            new["Substitute longest homopolymer"] = str(hp_ef.get("longest", 0))
            new["Substitute homopolymer count"] = str(hp_ef.get("count_ge2", 0))

        out_rows.append(new)
    return out_rows


def strand_rows_to_csv(rows: List[Dict[str, str]]) -> str:
    buf = io.StringIO()
    if not rows:
        return ""

    preferred = [
        "No.",
        "FBR",
        "Strand index",
        "Payload",
        "Filler",
        "RBR",
        "Full strand",
        "FBR length",
        "Index length",
        "Payload length",
        "Filler length",
        "RBR length",
        "Payload capacity",
        "Target total length",
        "Total length",
        "Original DNA length",
        "GC content",
        "Longest homopolymer",
        "Homopolymer count",
        "Substitute payload",
        "Substitute full strand",
        "Error scope",
        "Substitute payload length",
        "Substitute full length",
        "Substitute count",
        "Changed bases",
        "Error count",
        "Substitute GC content",
        "Substitute longest homopolymer",
        "Substitute homopolymer count",
    ]
    extra = []
    for row in rows:
        for key in row.keys():
            if key not in preferred and key not in extra:
                extra.append(key)
    fieldnames = [k for k in preferred if any(k in r for r in rows)] + extra

    writer = csv.DictWriter(buf, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    return buf.getvalue()


def strand_rows_to_preview(rows: List[Dict[str, str]], max_rows: int = 5) -> str:
    lines = []
    for row in rows[:max_rows]:
        full_strand = row.get("Full strand") or (
            row.get("FBR", "")
            + row.get("Strand index", "")
            + row.get("Payload", "")
            + row.get("Filler", "")
            + row.get("RBR", "")
        )
        lines.append(
            f"{row.get('No.', '')}. "
            f"FBR({len(row.get('FBR', ''))}) + "
            f"Index({row.get('Strand index', '')}) + "
            f"Payload({len(row.get('Payload', ''))}) + "
            f"Filler({len(row.get('Filler', ''))}) + "
            f"RBR({len(row.get('RBR', ''))}) = "
            f"Total({len(full_strand)})"
        )
        lines.append(full_strand)
        if row.get("Substitute full strand"):
            lines.append("Substitute full strand:")
            lines.append(row.get("Substitute full strand", ""))
        lines.append("")
    return "\n".join(lines)
