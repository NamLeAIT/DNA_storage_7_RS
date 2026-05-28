from __future__ import annotations

"""
Rigorous DNAStorageToolkit adapter for the compression-aware DNA storage app.

This module keeps the user's compression-first DNA encoder as the front end and
adds a downstream layer that follows the DNAStorageToolkit pipeline boundary:

    EncodedStrands.txt
      -> wetlab/noise simulation
      -> NoisyStrands.txt / UnderlyingClusters.txt
      -> clustering-like grouping
      -> ReconstructedStrands.txt
      -> recovered payload DNA for the app decoder

Two execution paths are supported:

1) Built-in rigorous mode
   Dependency-light implementation for Streamlit. It supports substitution,
   insertion, deletion, read dropout, index-aware clustering, and a fixed-length
   consensus reconstruction aligned to the prepared strand design.

2) External DNAStorageToolkit compatibility
   Export/import helpers and repo validation are included so that the real
   DNAStorageToolkit repository can be attached without changing the rest of the
   app. When the external repo is present, the app can at least run the official
   naive simulator and export/import the standard intermediate files. Clustering
   and trace reconstruction are kept as optional external shell steps because the
   official toolkit depends on compiled C++/make and pyspoa.
"""

import csv
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

BASES = "ACGT"


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------

def clean_dna(seq: str) -> str:
    return "".join(ch for ch in str(seq or "").upper() if ch in BASES)


def _safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except Exception:
        return default


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(str(value).strip())
    except Exception:
        return default


def _row_no(row: Dict[str, Any], fallback: int) -> str:
    return str(row.get("No.") or row.get("No") or row.get("strand_no") or fallback)


def _full_strand(row: Dict[str, Any]) -> str:
    for key in ("Full strand", "full_strand", "Full", "sequence", "Sequence"):
        seq = clean_dna(row.get(key, ""))
        if seq:
            return seq
    # fallback for rows that only store payload
    return clean_dna(row.get("Payload", ""))


def _payload(row: Dict[str, Any]) -> str:
    """Payload used as the reconstruction reference.

    For normal prepared strands this is the original Payload field.
    For Advanced Add Errors rows, the wet-lab input is already the error strand,
    so reconstruction accuracy should be measured against the payload segment of
    that error strand, not against the clean pre-error payload.
    """
    for key in ("Wet-lab reference payload", "Reconstruction reference payload", "Error payload", "Payload"):
        seq = clean_dna(row.get(key, ""))
        if seq:
            return seq
    return ""


def _is_toolkit_rs_direct(row: Dict[str, Any]) -> bool:
    """True when the row already is a DNAStorageToolkit 120-nt indexed strand.

    In this mode the full 120 nt must be passed to the Toolkit RS decoder.
    The first 8 nt are the Toolkit column index, but they are also part of
    the encoded Toolkit strand; they must not be trimmed away as an app index.
    """
    flag = str(row.get("Toolkit RS direct", row.get("toolkit_rs_direct", ""))).strip().lower()
    if flag in {"1", "true", "yes", "y"}:
        return True
    return bool(row.get("Toolkit column index")) and len(_full_strand(row)) == 120


def _hamming_distance(a: str, b: str) -> int:
    a = clean_dna(a)
    b = clean_dna(b)
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] != b[i]) + abs(len(a) - len(b))


def edit_distance(a: str, b: str, max_cutoff: Optional[int] = None) -> int:
    """Small dependency-free Levenshtein distance with optional early cutoff."""
    a = clean_dna(a)
    b = clean_dna(b)
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    if len(a) < len(b):
        a, b = b, a
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        current = [i]
        row_min = current[0]
        for j, cb in enumerate(b, 1):
            insert = current[j - 1] + 1
            delete = previous[j] + 1
            replace = previous[j - 1] + (ca != cb)
            val = min(insert, delete, replace)
            current.append(val)
            if val < row_min:
                row_min = val
        previous = current
        if max_cutoff is not None and row_min > max_cutoff:
            return max_cutoff + 1
    return previous[-1]


def _qgrams(seq: str, q: int = 5) -> set[str]:
    seq = clean_dna(seq)
    if len(seq) < q:
        return {seq} if seq else set()
    return {seq[i:i + q] for i in range(0, len(seq) - q + 1)}


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / max(1, len(a | b))


# ---------------------------------------------------------------------------
# Toolkit-compatible file formats
# ---------------------------------------------------------------------------

def export_encoded_strands(strand_rows: Sequence[Dict[str, Any]], path: str | Path) -> str:
    """Write DNAStorageToolkit-style EncodedStrands.txt: one full strand per line."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        for row in strand_rows:
            seq = _full_strand(row)
            if seq:
                f.write(seq + "\n")
    return str(p)


def write_underlying_clusters(reads: Sequence[Dict[str, Any]], path: str | Path) -> str:
    """Write Toolkit naive simulator cluster format: CLUSTER i then reads."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    grouped: Dict[str, List[str]] = defaultdict(list)
    for read in reads:
        src = str(read.get("Source No.", ""))
        seq = clean_dna(read.get("Read sequence", ""))
        if src and seq:
            grouped[src].append(seq)
    with p.open("w", encoding="utf-8") as f:
        for src in sorted(grouped, key=lambda x: _safe_int(x, 10**9)):
            f.write(f"CLUSTER {src}\n")
            for seq in grouped[src]:
                f.write(seq + "\n")
    return str(p)


def write_noisy_strands(reads: Sequence[Dict[str, Any]], path: str | Path, *, shuffle: bool = True, seed: int = 1) -> str:
    """Write Toolkit NoisyStrands.txt: shuffled reads without cluster labels."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    seqs = [clean_dna(r.get("Read sequence", "")) for r in reads if clean_dna(r.get("Read sequence", ""))]
    if shuffle:
        rng = random.Random(int(seed))
        rng.shuffle(seqs)
    with p.open("w", encoding="utf-8") as f:
        for seq in seqs:
            f.write(seq + "\n")
    return str(p)


def write_clustered_strands(clusters: Dict[str, Sequence[str]], path: str | Path) -> str:
    """
    Write a reconstruction-friendly clustered file.

    The official recon.py expects a first line to discard, then integer cluster
    delimiters followed by reads. This format is intentionally simple and is also
    easy for the app to parse.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        f.write("# ClusteredStrands generated by compression-aware adapter\n")
        for src in sorted(clusters, key=lambda x: _safe_int(x, 10**9)):
            f.write(str(src) + "\n")
            for seq in clusters[src]:
                seq = clean_dna(seq)
                if seq:
                    f.write(seq + "\n")
    return str(p)


def write_reconstructed_strands(rows: Sequence[Dict[str, Any]], path: str | Path) -> str:
    """Write Toolkit-style ReconstructedStrands.txt: one reconstructed full strand per line."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        for row in rows:
            seq = clean_dna(row.get("Reconstructed full strand", ""))
            if seq:
                f.write(seq + "\n")
    return str(p)


def rows_to_csv(rows: Iterable[Dict[str, Any]]) -> str:
    rows = list(rows)
    if not rows:
        return ""
    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Simulation: follows DNAStorageToolkit naive/noise.py semantics
# ---------------------------------------------------------------------------

def _mutate_naive_toolkit(
    seq: str,
    *,
    sub_rate: float,
    del_rate: float,
    ins_rate: float,
    rng: random.Random,
) -> Tuple[str, Dict[str, int], List[Dict[str, Any]]]:
    """
    Match the official naive/noise.py logic closely and keep a per-base event log.

    For each base w:
      if r < sub: append random base
      elif r < sub + ins: append random base then w
      elif r > sub + ins + del: append w
      else: deletion

    Note: the official simulator may substitute with the same base because it
    draws from A/G/T/C without excluding w; this implementation follows that.
    """
    seq = clean_dna(seq)
    out: List[str] = []
    events: List[Dict[str, Any]] = []
    subs = dels = inss = 0
    read_pos = 0
    for pos, w in enumerate(seq, start=1):
        r = rng.random()
        if r < sub_rate:
            new_base = rng.choice(BASES)
            out.append(new_base)
            read_pos += 1
            subs += 1
            events.append({
                "position_original": pos,
                "position_read": read_pos,
                "operation": "substitution",
                "from_base": w,
                "to_base": new_base,
            })
        elif r < sub_rate + ins_rate:
            ins_base = rng.choice(BASES)
            out.append(ins_base)
            read_pos += 1
            events.append({
                "position_original": pos,
                "position_read": read_pos,
                "operation": "insertion",
                "from_base": "",
                "to_base": ins_base,
            })
            out.append(w)
            read_pos += 1
            inss += 1
        elif r > sub_rate + ins_rate + del_rate:
            out.append(w)
            read_pos += 1
        else:
            dels += 1
            events.append({
                "position_original": pos,
                "position_read": read_pos + 1,
                "operation": "deletion",
                "from_base": w,
                "to_base": "",
            })
    return "".join(out), {
        "Substitution count": subs,
        "Insertion count": inss,
        "Deletion count": dels,
        "Error count": subs + inss + dels,
    }, events


def simulate_sequencing_reads(
    strand_rows: Sequence[Dict[str, Any]],
    *,
    coverage: int = 10,
    substitution_rate: float = 0.001,
    deletion_rate: float = 0.0,
    insertion_rate: float = 0.0,
    dropout_rate: float = 0.0,
    randomize_coverage: bool = False,
    seed: int = 1,
) -> Tuple[List[Dict[str, Any]], Dict[str, int | float]]:
    """Create noisy reads from prepared full strands using Toolkit-naive semantics."""
    rng = random.Random(int(seed))
    coverage = max(1, int(coverage))
    dropout_rate = max(0.0, min(1.0, float(dropout_rate)))
    sub_rate = max(0.0, min(1.0, float(substitution_rate)))
    del_rate = max(0.0, min(1.0, float(deletion_rate)))
    ins_rate = max(0.0, min(1.0, float(insertion_rate)))

    reads: List[Dict[str, Any]] = []
    total_sub = total_del = total_ins = 0
    total_template_nt = 0
    attempted_reads = 0
    dropped_reads = 0

    for row_idx, row in enumerate(strand_rows, start=1):
        source_no = _row_no(row, row_idx)
        full = _full_strand(row)
        if not full:
            continue

        n_reads = coverage
        if randomize_coverage:
            n_reads = max(1, int(round(rng.uniform(max(1, coverage * 0.5), max(1, coverage * 1.5)))))

        for copy_no in range(1, n_reads + 1):
            attempted_reads += 1
            if rng.random() < dropout_rate:
                dropped_reads += 1
                continue
            read_rng = random.Random(f"{seed}|{source_no}|{copy_no}|{sub_rate}|{del_rate}|{ins_rate}")
            seq, counts, events = _mutate_naive_toolkit(
                full,
                sub_rate=sub_rate,
                del_rate=del_rate,
                ins_rate=ins_rate,
                rng=read_rng,
            )
            total_sub += int(counts["Substitution count"])
            total_del += int(counts["Deletion count"])
            total_ins += int(counts["Insertion count"])
            total_template_nt += len(full)
            read_id = f"r{len(reads) + 1:06d}"
            for ev in events:
                ev["read_id"] = read_id
                ev["source_no"] = source_no
                ev["copy_no"] = str(copy_no)
            reads.append({
                "Read ID": read_id,
                "Source No.": source_no,
                "Copy No.": str(copy_no),
                "Read sequence": seq,
                "Original full strand": full,
                "Substitution count": str(counts["Substitution count"]),
                "Insertion count": str(counts["Insertion count"]),
                "Deletion count": str(counts["Deletion count"]),
                "Error count": str(counts["Error count"]),
                "Read length": str(len(seq)),
                "Template length": str(len(full)),
                "Error events": json.dumps(events, ensure_ascii=False),
                "Event preview": "; ".join(
                    f"{e['operation']}@{e['position_original']}" for e in events[:8]
                ) + ("; ..." if len(events) > 8 else ""),
            })

    source_with_reads = {str(r["Source No."]) for r in reads}
    metrics: Dict[str, int | float] = {
        "input_strands": len(strand_rows),
        "coverage": coverage,
        "attempted_reads": attempted_reads,
        "dropped_reads": dropped_reads,
        "reads_generated": len(reads),
        "strands_with_reads": len(source_with_reads),
        "lost_strands": max(0, len(strand_rows) - len(source_with_reads)),
        "total_substitutions": total_sub,
        "total_insertions": total_ins,
        "total_deletions": total_del,
        "total_errors": total_sub + total_ins + total_del,
        "observed_error_rate": ((total_sub + total_ins + total_del) / total_template_nt) if total_template_nt else 0.0,
    }
    return reads, metrics


# ---------------------------------------------------------------------------
# Clustering and reconstruction
# ---------------------------------------------------------------------------

def _prefix_len_from_template(row: Dict[str, Any]) -> int:
    # Toolkit RS strands already start with an 8-nt column index.  Use that
    # index for clustering, but do not remove it during payload extraction.
    if _is_toolkit_rs_direct(row):
        fbr = clean_dna(row.get("FBR", ""))
        idx = clean_dna(row.get("Index", row.get("Strand index", "")))
        if fbr or idx:
            return min(len(fbr) + len(idx), len(_full_strand(row)))
        return min(12, len(_full_strand(row)))
    fbr_len = _safe_int(row.get("FBR length"), len(clean_dna(row.get("FBR", ""))))
    idx_len = _safe_int(row.get("Index length"), len(clean_dna(row.get("Strand index", ""))))
    if fbr_len + idx_len > 0:
        return fbr_len + idx_len
    return min(28, len(_full_strand(row)))


def cluster_reads(
    strand_rows: Sequence[Dict[str, Any]],
    reads: Sequence[Dict[str, Any]],
    *,
    method: str = "index_aware",
    q: int = 5,
    max_prefix_edit: int = 6,
) -> Tuple[Dict[str, List[str]], Dict[str, int | float]]:
    """
    Group noisy reads into per-strand clusters.

    method:
      - oracle_source: uses simulator Source No.; best-case control only.
      - index_aware: assigns each read to the nearest designed prefix. This is
        practical for this app because prepared strands include designed indices.
      - qgram_nearest: assigns by whole-strand q-gram Jaccard similarity.
    """
    method = str(method or "index_aware")
    templates: List[Tuple[str, str, str, int, set[str]]] = []
    for i, row in enumerate(strand_rows, start=1):
        no = _row_no(row, i)
        full = _full_strand(row)
        if not full:
            continue
        plen = _prefix_len_from_template(row)
        templates.append((no, full, full[:plen], plen, _qgrams(full, q=q)))

    clusters: Dict[str, List[str]] = defaultdict(list)
    assigned = 0
    unassigned = 0
    ambiguous = 0

    if method == "oracle_source":
        valid = {t[0] for t in templates}
        for read in reads:
            src = str(read.get("Source No.", ""))
            seq = clean_dna(read.get("Read sequence", ""))
            if src in valid and seq:
                clusters[src].append(seq)
                assigned += 1
            elif seq:
                unassigned += 1
        return dict(clusters), {"cluster_method": 0, "assigned_reads": assigned, "unassigned_reads": unassigned, "ambiguous_reads": 0}

    for read in reads:
        seq = clean_dna(read.get("Read sequence", ""))
        if not seq:
            continue

        scored: List[Tuple[float, int, str]] = []
        if method == "qgram_nearest":
            rq = _qgrams(seq, q=q)
            for no, full, _prefix, _plen, tq in templates:
                jac = _jaccard(rq, tq)
                # lower score is better; edit distance is used as tie-breaker
                approx_ed = abs(len(seq) - len(full))
                scored.append((-jac, approx_ed, no))
        else:  # index_aware default
            for no, _full, prefix, plen, _tq in templates:
                # Compare the designed prefix at its exact length.
                #
                # Earlier versions used `plen + max_prefix_edit` bases from the
                # read to tolerate indels near the index. That made a clean read
                # score 6 rather than 0 because the extra payload bases had to be
                # deleted during edit-distance matching. Many different strands
                # could then tie, so even with noise=0 reads were assigned to the
                # wrong cluster and reconstructed payload accuracy dropped below
                # 1.0. Exact-length prefix matching is the correct default for
                # prepared strands; indel tolerance should come from later read
                # alignment, not from adding payload bases to the clustering key.
                read_prefix = seq[: max(1, min(len(seq), plen))]
                d = edit_distance(read_prefix, prefix, max_cutoff=max_prefix_edit + 4)
                scored.append((float(d), abs(len(read_prefix) - len(prefix)), no))

        if not scored:
            unassigned += 1
            continue
        scored.sort(key=lambda x: (x[0], x[1], _safe_int(x[2], 10**9)))
        best = scored[0]
        second = scored[1] if len(scored) > 1 else None
        if method == "index_aware" and best[0] > max_prefix_edit + 4:
            unassigned += 1
            continue
        if second is not None and best[0] == second[0] and best[1] == second[1]:
            ambiguous += 1
        clusters[best[2]].append(seq)
        assigned += 1

    return dict(clusters), {
        "cluster_method": {"oracle_source": 0, "index_aware": 1, "qgram_nearest": 2}.get(method, 1),
        "assigned_reads": assigned,
        "unassigned_reads": unassigned,
        "ambiguous_reads": ambiguous,
    }


def _align_to_reference(read: str, ref: str) -> Tuple[List[Optional[str]], Dict[int, List[str]]]:
    """
    Global align read to reference. Return base calls per reference position and
    insertion calls before reference positions.

    The reconstruction uses reference-position votes and ignores insertions by
    default because the expected full-strand length is known from design.
    """
    read = clean_dna(read)
    ref = clean_dna(ref)
    n, m = len(read), len(ref)
    # DP cost matrix; strings are short (~120 nt), so this is fine for app scale.
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    bt = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = i
        bt[i][0] = "U"  # insertion relative to ref
    for j in range(1, m + 1):
        dp[0][j] = j
        bt[0][j] = "L"  # deletion in read
    for i in range(1, n + 1):
        ri = read[i - 1]
        for j in range(1, m + 1):
            rj = ref[j - 1]
            diag = dp[i - 1][j - 1] + (ri != rj)
            up = dp[i - 1][j] + 1
            left = dp[i][j - 1] + 1
            best = min(diag, up, left)
            dp[i][j] = best
            if best == diag:
                bt[i][j] = "D"
            elif best == left:
                bt[i][j] = "L"
            else:
                bt[i][j] = "U"

    calls: List[Optional[str]] = [None] * m
    insertions: Dict[int, List[str]] = defaultdict(list)
    i, j = n, m
    while i > 0 or j > 0:
        step = bt[i][j]
        if i > 0 and j > 0 and step == "D":
            calls[j - 1] = read[i - 1]
            i -= 1
            j -= 1
        elif j > 0 and (i == 0 or step == "L"):
            calls[j - 1] = None
            j -= 1
        else:
            # read base inserted before ref position j
            if i > 0:
                insertions[j].append(read[i - 1])
            i -= 1
    return calls, insertions


def _majority_base(counts: Counter[str], fallback: str = "A") -> str:
    if not counts:
        return fallback if fallback in BASES else "A"
    return max(BASES, key=lambda b: (counts.get(b, 0), -BASES.index(b)))


def _consensus_against_template(reads: Sequence[str], template_full: str) -> str:
    """Consensus at known template length using global alignment of reads."""
    ref = clean_dna(template_full)
    if not reads:
        return ""
    if not ref:
        ref = max((clean_dna(r) for r in reads), key=len, default="")
    votes: List[Counter[str]] = [Counter() for _ in range(len(ref))]
    for read in reads:
        seq = clean_dna(read)
        if not seq:
            continue
        # Exact length reads are common under substitution-only simulation.
        if len(seq) == len(ref):
            for i, b in enumerate(seq):
                if b in BASES:
                    votes[i][b] += 1
            continue
        calls, _insertions = _align_to_reference(seq, ref)
        for i, base in enumerate(calls):
            if base is not None and base in BASES:
                votes[i][base] += 1
    return "".join(_majority_base(votes[i], fallback=ref[i] if i < len(ref) else "A") for i in range(len(ref)))


def _payload_from_full_strand(full: str, template: Dict[str, Any]) -> str:
    full = clean_dna(full)
    if _is_toolkit_rs_direct(template):
        # For Reed-Solomon direct strands, recover only the coding core:
        # Index + Payload + RS parity. FBR/RBR are sequencing wrappers.
        fbr = clean_dna(template.get("FBR", ""))
        rbr = clean_dna(template.get("RBR", ""))
        if fbr and rbr and len(full) >= len(fbr) + len(rbr):
            return full[len(fbr):len(full) - len(rbr)]
        return clean_dna(
            template.get("Index", template.get("Strand index", ""))
            + template.get("Payload", "")
            + template.get("RS parity", "")
        ) or full
    fbr_len = _safe_int(template.get("FBR length"), len(clean_dna(template.get("FBR", ""))))
    idx_len = _safe_int(template.get("Index length"), len(clean_dna(template.get("Strand index", ""))))
    payload_len = _safe_int(template.get("Payload length"), len(clean_dna(template.get("Payload", ""))))
    start = max(0, fbr_len + idx_len)
    end = min(len(full), start + max(0, payload_len))
    return full[start:end]


def reconstruct_consensus_from_reads(
    strand_rows: Sequence[Dict[str, Any]],
    reads: Sequence[Dict[str, Any]],
    *,
    cluster_method: str = "index_aware",
    q: int = 5,
    max_prefix_edit: int = 6,
) -> Tuple[List[Dict[str, Any]], str, Dict[str, int | float]]:
    """Cluster reads, reconstruct full strands, and join recovered payload DNA."""
    clusters, cluster_metrics = cluster_reads(
        strand_rows,
        reads,
        method=cluster_method,
        q=q,
        max_prefix_edit=max_prefix_edit,
    )

    out_rows: List[Dict[str, Any]] = []
    payload_parts: List[str] = []
    full_mismatches = 0
    payload_mismatches = 0
    total_full_nt = 0
    total_payload_nt = 0

    for row_idx, template in enumerate(strand_rows, start=1):
        source_no = _row_no(template, row_idx)
        source_reads = clusters.get(source_no, [])
        if not source_reads:
            continue

        original_full = _full_strand(template)
        original_payload = _payload(template)
        reconstructed_full = _consensus_against_template(source_reads, original_full)
        reconstructed_payload = _payload_from_full_strand(reconstructed_full, template)

        fd = _hamming_distance(original_full, reconstructed_full)
        pd = _hamming_distance(original_payload, reconstructed_payload)
        full_mismatches += fd
        payload_mismatches += pd
        total_full_nt += max(len(original_full), len(reconstructed_full))
        total_payload_nt += max(len(original_payload), len(reconstructed_payload))
        payload_parts.append(reconstructed_payload)

        out = dict(template)
        out.update({
            "Reads used": str(len(source_reads)),
            "Reconstructed full strand": reconstructed_full,
            "Reconstructed payload": reconstructed_payload,
            "Consensus full mismatches": str(fd),
            "Consensus payload mismatches": str(pd),
        })
        out_rows.append(out)

    metrics: Dict[str, int | float] = {
        **cluster_metrics,
        "reconstructed_strands": len(out_rows),
        "missing_strands": max(0, len(strand_rows) - len(out_rows)),
        "full_mismatches": full_mismatches,
        "payload_mismatches": payload_mismatches,
        "full_consensus_accuracy": 1.0 - (full_mismatches / total_full_nt) if total_full_nt else 0.0,
        "payload_consensus_accuracy": 1.0 - (payload_mismatches / total_payload_nt) if total_payload_nt else 0.0,
    }
    return out_rows, "".join(payload_parts), metrics


# ---------------------------------------------------------------------------
# External DNAStorageToolkit hooks
# ---------------------------------------------------------------------------

@dataclass
class ExternalToolkitStatus:
    ok: bool
    repo_path: str
    message: str
    found: Dict[str, bool]


def validate_external_toolkit_repo(repo_path: str | Path) -> ExternalToolkitStatus:
    root = Path(str(repo_path)).expanduser()
    required = {
        "noise.py": root / "2-simulating_wetlab" / "naive" / "noise.py",
        "shuffle.py": root / "2-simulating_wetlab" / "naive" / "shuffle.py",
        "recon.py": root / "4-reconstruction" / "recon.py",
        "clustering_dir": root / "3-clustering",
        "codec.py": root / "1-encoding-decoding" / "codec.py",
    }
    found = {name: path.exists() for name, path in required.items()}
    ok = bool(root.exists() and found["noise.py"] and found["shuffle.py"])
    missing = [name for name, exists in found.items() if not exists]
    if ok:
        msg = "External DNAStorageToolkit repo detected. Naive wetlab simulator is available; clustering/reconstruction depend on its local dependencies."
    else:
        msg = "Not a valid DNAStorageToolkit path" + (f"; missing: {', '.join(missing)}" if missing else ".")
    return ExternalToolkitStatus(ok=ok, repo_path=str(root), message=msg, found=found)


def run_external_naive_simulator(
    repo_path: str | Path,
    strand_rows: Sequence[Dict[str, Any]],
    work_dir: str | Path,
    *,
    coverage: int = 10,
    substitution_rate: float = 0.01,
    deletion_rate: float = 0.01,
    insertion_rate: float = 0.01,
) -> Dict[str, Any]:
    """
    Run official DNAStorageToolkit naive/noise.py + shuffle.py.

    This does not run compiled clustering or pyspoa reconstruction automatically;
    it prepares the standard files so the user can run or inspect them. The app's
    built-in reconstruction can still use the imported reads afterwards.
    """
    status = validate_external_toolkit_repo(repo_path)
    if not status.ok:
        raise RuntimeError(status.message)

    root = Path(status.repo_path)
    wd = Path(work_dir)
    wd.mkdir(parents=True, exist_ok=True)
    encoded = wd / "EncodedStrands.txt"
    underlying = wd / "UnderlyingClusters.txt"
    noisy = wd / "NoisyStrands.txt"
    export_encoded_strands(strand_rows, encoded)

    noise_py = root / "2-simulating_wetlab" / "naive" / "noise.py"
    shuffle_py = root / "2-simulating_wetlab" / "naive" / "shuffle.py"
    cmd_noise = [
        sys.executable,
        str(noise_py),
        "--N", str(int(coverage)),
        "--subs", str(float(substitution_rate)),
        "--dels", str(float(deletion_rate)),
        "--inss", str(float(insertion_rate)),
        "--i", str(encoded),
        "--o", str(underlying),
    ]
    cmd_shuffle = [sys.executable, str(shuffle_py), str(underlying), str(noisy)]
    p1 = subprocess.run(cmd_noise, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p1.returncode != 0:
        raise RuntimeError("DNAStorageToolkit noise.py failed: " + (p1.stderr or p1.stdout)[-2000:])
    p2 = subprocess.run(cmd_shuffle, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p2.returncode != 0:
        raise RuntimeError("DNAStorageToolkit shuffle.py failed: " + (p2.stderr or p2.stdout)[-2000:])
    return {
        "repo_status": status.message,
        "work_dir": str(wd),
        "encoded_path": str(encoded),
        "underlying_path": str(underlying),
        "noisy_path": str(noisy),
        "noise_stdout": p1.stdout,
        "shuffle_stdout": p2.stdout,
    }


def read_noisy_strands(path: str | Path) -> List[str]:
    p = Path(path)
    if not p.exists():
        return []
    out: List[str] = []
    for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        seq = clean_dna(line)
        if seq:
            out.append(seq)
    return out


def wrap_noisy_sequences_as_reads(seqs: Sequence[str]) -> List[Dict[str, Any]]:
    return [{"Read ID": f"ext{i+1:06d}", "Read sequence": clean_dna(seq)} for i, seq in enumerate(seqs) if clean_dna(seq)]


# ---------------------------------------------------------------------------
# One-call pipeline used by Streamlit
# ---------------------------------------------------------------------------

def run_builtin_toolkit_pipeline(
    strand_rows: Sequence[Dict[str, Any]],
    *,
    coverage: int = 10,
    substitution_rate: float = 0.001,
    deletion_rate: float = 0.0,
    insertion_rate: float = 0.0,
    dropout_rate: float = 0.0,
    randomize_coverage: bool = False,
    cluster_method: str = "index_aware",
    q: int = 5,
    max_prefix_edit: int = 6,
    seed: int = 1,
    export_dir: str | Path | None = None,
) -> Dict[str, Any]:
    reads, sim_metrics = simulate_sequencing_reads(
        strand_rows,
        coverage=coverage,
        substitution_rate=substitution_rate,
        deletion_rate=deletion_rate,
        insertion_rate=insertion_rate,
        dropout_rate=dropout_rate,
        randomize_coverage=randomize_coverage,
        seed=seed,
    )
    recon_rows, payload_dna, rec_metrics = reconstruct_consensus_from_reads(
        strand_rows,
        reads,
        cluster_method=cluster_method,
        q=q,
        max_prefix_edit=max_prefix_edit,
    )
    metrics = {**sim_metrics, **rec_metrics}

    paths: Dict[str, str] = {}
    if export_dir is not None:
        out = Path(export_dir)
        out.mkdir(parents=True, exist_ok=True)
        paths["EncodedStrands.txt"] = export_encoded_strands(strand_rows, out / "EncodedStrands.txt")
        paths["UnderlyingClusters.txt"] = write_underlying_clusters(reads, out / "UnderlyingClusters.txt")
        paths["NoisyStrands.txt"] = write_noisy_strands(reads, out / "NoisyStrands.txt", seed=seed)
        clusters, _cluster_metrics = cluster_reads(strand_rows, reads, method=cluster_method, q=q, max_prefix_edit=max_prefix_edit)
        paths["ClusteredStrands.txt"] = write_clustered_strands(clusters, out / "ClusteredStrands.txt")
        paths["ReconstructedStrands.txt"] = write_reconstructed_strands(recon_rows, out / "ReconstructedStrands.txt")

    return {
        "reads": reads,
        "reconstructed_rows": recon_rows,
        "payload_dna": payload_dna,
        "metrics": metrics,
        "paths": paths,
        "reads_csv": rows_to_csv(reads),
        "reconstructed_csv": rows_to_csv(recon_rows),
    }
