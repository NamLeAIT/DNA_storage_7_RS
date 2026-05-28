from __future__ import annotations

import gzip
import io
import lzma
import bz2
import zipfile
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

try:
    from PIL import Image
except Exception:
    Image = None

import dna_codec
import new_design_codec
import toolkit_rs_codec
import rs_concatenated_codec
from utils_bits_v2 import bytes_to_bitstring, bitstring_to_bytes, detect_magic
from config import MAPPING_OPTIONS, IMAGE_KINDS


def mapping_to_config(mapping_name: str) -> Dict[str, Any]:
    if mapping_name == "Simple Mapping":
        return {
            "mode": "SIMPLE",
            "scheme_name": "RINF_B16",
            "init_dimer": "TA",
            "whiten": False,
        }
    return {
        "mode": "TABLE",
        "scheme_name": mapping_name,
        "init_dimer": "TA",
        "whiten": False,
    }

def encode_bytes_to_dna(data: bytes, mapping_name: str) -> Tuple[str, str, Dict[str, Any]]:
    bits = bytes_to_bitstring(data)
    if bits == "":
        bits = "0"

    if mapping_name == "Reed-Solomon":
        result = rs_concatenated_codec.encode(data)
        meta = dict(result.meta)
        meta.update({
            "mapping": mapping_name,
            "mode": "REED_SOLOMON",
            "bits_len": len(bits),
            "bytes_len": len(data),
        })
        return result.dna, result.bits, meta

    if mapping_name == "New Design":
        result = new_design_codec.encode_new_design(data)
        meta = dict(result.get("meta", {}))
        meta.update({
            "mapping": mapping_name,
            "mode": "NEW_DESIGN",
            "bits_len": len(bits),
            "bytes_len": len(data),
            "blocks": result.get("blocks", []),
        })
        return result["dna"], bits, meta

    if mapping_name == "Toolkit RS Baseline":
        result = toolkit_rs_codec.encode_toolkit_rs(data, redundancy_pct=15.0)
        meta = dict(result.meta)
        meta.update({
            "bits_len": len(bits),
            "bytes_len": len(data),
            "toolkit_strands": result.strands,
        })
        return result.dna, bits, meta

    cfg = mapping_to_config(mapping_name)
    dna, digits = dna_codec.encode_bits_to_dna(
        bits,
        scheme_name=cfg["scheme_name"],
        mode=cfg["mode"],
        seed="rn",
        init_dimer=cfg["init_dimer"],
        prepend_one=True,
        whiten=cfg["whiten"],
        target_gc=0.50,
        w_gc=0.0,
        w_motif=0.0,
        ks=(4, 6),
    )
    meta = {
        "mapping": mapping_name,
        "mode": cfg["mode"],
        "scheme_name": cfg["scheme_name"],
        "init_dimer": "TA",
        "bits_len": len(bits),
        "digits_len": len(digits) if isinstance(digits, list) else None,
        "bytes_len": len(data),
    }
    return dna, bits, meta

def decode_dna_with_mapping(dna: str, mapping_name: str, codec_meta: Optional[Dict[str, Any]] = None) -> Tuple[bytes, str, Dict[str, Any]]:
    if mapping_name == "Toolkit RS Baseline":
        codec_meta = codec_meta or {}
        data, meta = toolkit_rs_codec.decode_toolkit_rs(
            dna,
            file_size=codec_meta.get("file_size") or codec_meta.get("bytes_len"),
            data_columns=codec_meta.get("data_columns"),
            parity_columns=codec_meta.get("parity_columns"),
            redundancy_pct=float(codec_meta.get("redundancy_pct_final", 15.0)),
        )
        bits = bytes_to_bitstring(data)
        meta.update({"bits_len": len(bits), "bytes_len": len(data)})
        return data, bits, meta

    if mapping_name == "Reed-Solomon":
        data, meta = rs_concatenated_codec.decode(dna, codec_meta or {})
        bits = bytes_to_bitstring(data)
        meta.update({"bits_len": len(bits), "bytes_len": len(data), "mapping": mapping_name})
        return data, bits, meta

    if mapping_name == "New Design":
        result = new_design_codec.decode_new_design(dna)

        ranked_candidates = []
        for cand in result.get("candidates", []):
            cdata = cand.get("data", b"")
            magic = detect_magic(cdata)
            score = -float(cand.get("distance", 0))
            note = "No recognizable file signature"
            valid = False
            if magic:
                score += 100.0 * float(magic.confidence)
                valid, note = validate_container(cdata, magic.kind)
                if valid:
                    score += 50.0
                else:
                    score -= 10.0
            ranked_candidates.append({
                "candidate": cand,
                "score": score,
                "magic": magic,
                "valid": valid,
                "note": note,
            })

        if ranked_candidates:
            ranked_candidates.sort(key=lambda x: x["score"], reverse=True)
            chosen = ranked_candidates[0]
            data = chosen["candidate"].get("data", b"")
            chosen_rank = int(chosen["candidate"].get("rank", 1))
        else:
            data = result["best_data"]
            chosen = None
            chosen_rank = 1

        bits = bytes_to_bitstring(data)
        meta = dict(result.get("meta", {}))
        candidate_outputs = []
        for item in ranked_candidates:
            cdata = item["candidate"].get("data", b"")
            cmagic = item["magic"]
            candidate_outputs.append({
                **{k: v for k, v in item["candidate"].items() if k != "data"},
                "score": item["score"],
                "detected": cmagic.kind if cmagic else "unknown",
                "ext": cmagic.ext if cmagic else ".bin",
                "valid_file": bool(item["valid"]),
                "note": item["note"],
                "data": cdata,
            })

        meta.update({
            "mapping": mapping_name,
            "bits_len": len(bits),
            "bytes_len": len(data),
            "repair_status": result.get("status"),
            "corrected_bases": result.get("corrected_bases", 0),
            "candidate_count": result.get("candidate_count", 0),
            "chosen_candidate": chosen_rank,
            "candidates": [
                {k: v for k, v in cand.items() if k != "data"}
                for cand in candidate_outputs
            ],
            "candidate_outputs": candidate_outputs,
        })
        return data, bits, meta

    cfg = mapping_to_config(mapping_name)
    bits, digits = dna_codec.decode_dna_to_bits(
        dna,
        scheme_name=cfg["scheme_name"],
        mode=cfg["mode"],
        seed="rn",
        init_dimer=cfg["init_dimer"],
        remove_leading_one=True,
        whiten=cfg["whiten"],
        target_gc=0.50,
        w_gc=0.0,
        w_motif=0.0,
        ks=(4, 6),
    )
    data, pad_bits = bitstring_to_bytes(bits, pad_to_byte=True)
    meta = {
        "mapping": mapping_name,
        "bits_len": len(bits),
        "bytes_len": len(data),
        "pad_bits_to_byte": pad_bits,
        "digits_len": len(digits) if isinstance(digits, list) else None,
    }
    return data, bits, meta

def validate_container(data: bytes, magic_kind: str) -> Tuple[bool, str]:
    """Lightweight validation beyond magic signature."""
    try:
        if magic_kind == "zip" or magic_kind in {"docx", "pptx", "xlsx", "epub"}:
            with zipfile.ZipFile(io.BytesIO(data), "r") as zf:
                bad = zf.testzip()
                if bad is not None:
                    return False, f"ZIP test failed at {bad}"
            return True, "ZIP container opened successfully"
        if magic_kind == "gzip":
            gzip.decompress(data)
            return True, "GZIP decompressed successfully"
        if magic_kind == "xz":
            lzma.decompress(data, format=lzma.FORMAT_XZ)
            return True, "XZ decompressed successfully"
        if magic_kind == "bz2":
            bz2.decompress(data)
            return True, "BZ2 decompressed successfully"
        if magic_kind in IMAGE_KINDS and Image is not None:
            img = Image.open(io.BytesIO(data))
            img.verify()
            return True, "Image verified successfully"
        return True, "Magic signature accepted"
    except Exception as e:
        return False, str(e)

def blind_decode_dna(dna_text: str) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """
    Decode DNA by trying available mappings.

    Important guard: New Design decoding can be intentionally expensive because it
    enumerates block candidates for repair. If a non-New-Design mapping already
    produces a verified self-describing file/container, we stop before trying New
    Design. This avoids treating ordinary mapped DNA as a repair-coded New Design
    stream and prevents long/hanging auto-detection.
    """
    dna = dna_codec.clean_dna_text(dna_text)
    rows: List[Dict[str, Any]] = []
    best: Optional[Dict[str, Any]] = None

    for mapping in MAPPING_OPTIONS:
        if mapping == "Toolkit RS Baseline":
            rows.append({
                "Mapping": mapping,
                "Status": "Manual only",
                "Magic": "—",
                "Ext": "—",
                "Confidence": 0.0,
                "Bytes": 0,
                "Score": 0.0,
                "Note": "Toolkit RS Baseline needs codec metadata such as original file size and data/parity columns; use Manual mapping after encoding in the same session.",
            })
            continue

        # Fast-stop before the expensive New Design candidate search if a normal
        # mapping has already produced a valid file/container.
        if mapping == "New Design" and best is not None and best.get("row", {}).get("Status") == "Valid":
            rows.append({
                "Mapping": mapping,
                "Status": "Skipped",
                "Magic": "—",
                "Ext": "—",
                "Confidence": 0.0,
                "Bytes": 0,
                "Score": 0.0,
                "Note": "Skipped because an earlier mapping already produced a valid self-describing file.",
            })
            break

        row: Dict[str, Any] = {"Mapping": mapping}
        try:
            data, bits, meta = decode_dna_with_mapping(dna, mapping)
            m = detect_magic(data)
            score = 0.0
            valid = False
            if data:
                score += 1.0
            if m:
                score += 10.0 * float(m.confidence)
                ok, note = validate_container(data, m.kind)
                valid = bool(ok)
                if ok:
                    score += 5.0
                else:
                    score -= 2.0
                row.update({
                    "Status": "Valid" if ok else "Weak",
                    "Magic": m.kind,
                    "Ext": m.ext,
                    "Confidence": m.confidence,
                    "Bytes": len(data),
                    "Score": score,
                    "Note": note,
                })
            else:
                row.update({
                    "Status": "No magic",
                    "Magic": "—",
                    "Ext": "—",
                    "Confidence": 0.0,
                    "Bytes": len(data),
                    "Score": score,
                    "Note": "No recognizable file/container signature",
                })

            candidate = {
                "mapping": mapping,
                "data": data,
                "bits": bits,
                "meta": meta,
                "magic": m,
                "score": score,
                "row": row,
            }
            if m is not None and (best is None or candidate["score"] > best["score"]):
                best = candidate


        except Exception as e:
            row.update({
                "Status": "Failed",
                "Magic": "—",
                "Ext": "—",
                "Confidence": 0.0,
                "Bytes": 0,
                "Score": -1.0,
                "Note": str(e)[:160],
            })
        rows.append(row)

    df = pd.DataFrame(rows)
    if best is None:
        raise ValueError("Auto-detection failed: no mapping produced a recognizable self-describing byte stream.")
    return best, df
