"""Pack the folded models of long proteins into one compact bundle per protein for the Giants tab.

MODELS_TSV lists one model per row (tab-separated, with a header line):
- protein: record id in COHORT_FASTA;
- scheme: the window scheme the model belongs to, w<size>_s<step> such as w2000_s1000 (the page's model-size buttons);
- start: 0-based position of the model's first residue in the cohort sequence;
- path: the model, a single-chain mmCIF with pLDDT in the B-factor column;
- optional mode (how the model was folded, shown under it in the page), ptm and ranking (its pTM and ranking
  score), and offsets: comma-separated shifts to try where the folded sequence and the cohort one differ by
  insertions (default 0).

Each model's residues are checked against the cohort sequence at start + offset (X / UNK match anything; up to
MAX_SUBSTITUTIONS mismatches) and encoded with backbone_codec (5 bytes per residue). A protein without models gets
an empty bundle, and the page shows Talea's calls and attention for it.

Writes OUT_DIR/bundles/<nn>.bin per protein in FASTA order: zlib(uint32 header length | header JSON | records), the
header holding the protein id, sequence and windows [scheme, start, end (0-based, half-open, cohort frame), record
offset, residues, mean pLDDT, pTM, ranking score, fold mode, substitutions]; OUT_DIR/bundles.json mapping ids to
files; and OUT_DIR/failures.json listing every model that did not read, match or encode.

usage: python3 encode_structures.py COHORT_FASTA MODELS_TSV OUT_DIR   (workers: SLURM_CPUS_ON_NODE, else all CPUs)
"""
import csv
import json
import os
import re
import struct
import sys
import zlib
from multiprocessing import Pool
from pathlib import Path

import backbone_codec

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
MAX_SUBSTITUTIONS = 3


def read_fasta(path):
    records, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name:
            records[name].append(line)
    return {name: "".join(parts) for name, parts in records.items()}


def number(value):
    return None if value in (None, "") else float(value)


def scheme_order(scheme):
    """Schemes named w<size>_s<step> by window size, then step; any other name after them."""
    match = re.fullmatch(r"w(\d+)_s(\d+)", scheme)
    return (0, int(match[1]), int(match[2]), "") if match else (1, 0, 0, scheme)


def encode_window(task):
    window, sequence, offsets = task
    try:
        model = backbone_codec.read_backbone(window["cif"])
    except (OSError, ValueError) as error:
        return window, None, f"read: {error}"
    residues = "".join(THREE_TO_ONE.get(name, "X") for name in model["resn"])
    best = None
    for offset in offsets:
        start = window["start"] + offset
        target = sequence[start:start + len(residues)]
        if start < 0 or len(target) != len(residues):
            continue
        substitutions = sum(1 for a, b in zip(residues, target) if a != b and "X" not in (a, b))
        if best is None or substitutions < best[1]:
            best = (start, substitutions)
    if best is None or best[1] > MAX_SUBSTITUTIONS:
        return window, None, f"sequence: {best}"
    try:
        blob = backbone_codec.encode(model)
    except ValueError as error:
        return window, None, f"encode: {error}"
    plddt = float(model["plddt"].mean())
    return {**window, "start": best[0], "end": best[0] + len(residues), "substitutions": best[1],
            "plddt": plddt}, blob, None


def main():
    cohort = read_fasta(sys.argv[1])
    out = Path(sys.argv[3])
    (out / "bundles").mkdir(parents=True, exist_ok=True)
    tasks = []
    for row in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
        window = {"protein": row["protein"], "scheme": row["scheme"], "start": int(row["start"]), "cif": row["path"],
                  "mode": row.get("mode") or "", "ptm": number(row.get("ptm")), "ranking": number(row.get("ranking"))}
        offsets = [int(value) for value in (row.get("offsets") or "0").split(",") if value]
        tasks.append((window, cohort[row["protein"]], offsets))
    print(f"{len(tasks)} models for {len(cohort)} proteins", flush=True)
    encoded, failures = {}, []
    with Pool(int(os.environ.get("SLURM_CPUS_ON_NODE") or os.cpu_count() or 1)) as pool:
        for done, (window, blob, error) in enumerate(pool.imap_unordered(encode_window, tasks, chunksize=16), 1):
            if error:
                failures.append({"protein": window["protein"], "scheme": window["scheme"], "start": window["start"],
                                 "cif": window["cif"], "error": error})
            else:
                encoded.setdefault(window["protein"], []).append((window, blob))
            if done % 2000 == 0:
                print(f"{done}/{len(tasks)}", flush=True)
    manifest = {}
    for index, protein in enumerate(cohort):
        items = sorted(encoded.get(protein, []), key=lambda item: (scheme_order(item[0]["scheme"]), item[0]["start"]))
        windows, records, offset = [], [], 0
        for window, blob in items:
            windows.append([window["scheme"], window["start"], window["end"], offset, window["end"] - window["start"],
                            round(window["plddt"], 1), window["ptm"] and round(window["ptm"], 3),
                            window["ranking"] and round(window["ranking"], 3), window["mode"],
                            window["substitutions"]])
            records.append(blob)
            offset += len(blob)
        header = json.dumps({"id": protein, "sequence": cohort[protein], "windows": windows},
                            separators=(",", ":")).encode()
        payload = struct.pack("<I", len(header)) + header + b"".join(records)
        name = f"{index:02d}.bin"
        (out / "bundles" / name).write_bytes(zlib.compress(payload, 9))
        manifest[protein] = {"file": name, "models": len(windows)}
    (out / "bundles.json").write_text(json.dumps(manifest, indent=1))
    (out / "failures.json").write_text(json.dumps(failures, indent=1))
    print(f"{sum(len(v) for v in encoded.values())} models encoded, {len(failures)} failed", flush=True)


if __name__ == "__main__":
    main()
