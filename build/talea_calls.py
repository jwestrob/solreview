"""Talea's region calls for the review bench: every repeat region with its units, positions 1-based inclusive.

Runs the Talea tree on PYTHONPATH (MPS) over every model of the named review sets, writes the full rows to ROWS_JSONL
(kept off the page) and, into the bench data, talea_calls.json with only where Talea places repeats: region
extents, period, Talea's region source and units, as sequence positions; a model numbered other than 1..L also gets
`resi`, the residue number at each position. Scores stay out of the page.

Rows that carry Talea's five-class type head also give each region `type` (from `region_type_probability`, the head
averaged over the region) and each call `protein_type` (from `predicted_teacher_mean_class_probability`, the head
averaged over the protein).

With --from-rows, ROWS_JSONL is read rather than written: rows from an earlier Talea run, matched to the models by
sequence SHA-256, so Talea does not run again. --talea-version then names the Talea that wrote them.

Calls for the other sets in the manifest are kept.

usage: PYTHONPATH=<Talea tree> python3 talea_calls.py BENCH_DIR ROWS_JSONL SET_ID=MODEL_DIR [SET_ID=MODEL_DIR ...]
       python3 talea_calls.py --from-rows --talea-version V BENCH_DIR ROWS_JSONL SET_ID=MODEL_DIR [...]
"""
import argparse
import hashlib
import json
from pathlib import Path

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def model_residues(path):
    numbers, sequence = [], []
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            numbers.append(int(line[22:26]))
            sequence.append(THREE_TO_ONE[line[17:20].strip()])
    return numbers, "".join(sequence)


def one_based(span):
    return {"start": int(span["start"]) + 1, "end": int(span["end"])}


def probabilities(values):
    return {name: round(float(value), 3) for name, value in values.items()}


def call_from_row(row):
    """The page's call for one Talea row: regions 1-based inclusive, plus the type head where the row has it."""
    regions = []
    for region in row.get("repeat_regions") or []:
        call = {
            **one_based(region),
            "period": round(float(region["period"]), 1),
            "source": region["source"],
            "units": [one_based(unit) for unit in region.get("units") or []],
        }
        if region.get("region_type_probability"):
            call["type"] = probabilities(region["region_type_probability"])
        regions.append(call)
    call = {"regions": regions}
    if row.get("talea_discovery_score") is not None:
        # Talea's score and its two halves: the residue model's z and the whole-protein model's z.
        call["scores"] = {
            "talea": round(float(row["talea_discovery_score"]), 3),
            "residue_z": round(float(row["teacher_state_z"]), 3),
            "whole_protein_z": round(float(row["global_z"]), 3),
        }
    if row.get("predicted_teacher_mean_class_probability"):
        call["protein_type"] = probabilities(row["predicted_teacher_mean_class_probability"])
    return call


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("bench", type=Path)
    parser.add_argument("rows", type=Path)
    parser.add_argument("sets", nargs="+", metavar="SET_ID=MODEL_DIR")
    parser.add_argument("--from-rows", action="store_true", help="read Talea rows from ROWS_JSONL instead of running Talea")
    parser.add_argument("--talea-version", help="with --from-rows: the Talea version that wrote the rows")
    args = parser.parse_args()
    if args.from_rows and not args.talea_version:
        parser.error("--from-rows needs --talea-version")

    model_dirs = dict(argument.split("=", 1) for argument in args.sets)
    manifest = json.loads((args.bench / "data" / "manifest.json").read_text())
    missing = sorted(set(model_dirs) - set(manifest["sets"]))
    if missing:
        raise SystemExit(f"not in the manifest: {', '.join(missing)}")
    calls_path = args.bench / "data" / "talea_calls.json"
    previous = json.loads(calls_path.read_text()) if calls_path.exists() else {"calls": {}}

    records, residue_numbers = [], {}
    for set_id in model_dirs:
        for item in manifest["sets"][set_id]:
            numbers, sequence = model_residues(Path(model_dirs[set_id]) / item["file"])
            # Talea positions index the modelled sequence; a model numbered other than 1..L also gets the
            # residue number of every position, so the page can show regions in the model's own numbering.
            if numbers != list(range(1, len(sequence) + 1)):
                residue_numbers[item["id"]] = numbers
            records.append((item["id"], sequence, hashlib.sha256(sequence.encode()).hexdigest()))

    calls = {}
    if args.from_rows:
        version = args.talea_version
        by_hash = {}
        for line in open(args.rows):
            row = json.loads(line)
            by_hash.setdefault(row["sequence_sha256"], row)
        unmatched = [record_id for record_id, _, digest in records if digest not in by_hash]
        if unmatched:
            raise SystemExit(f"{len(unmatched)} models have no row in {args.rows}, first {unmatched[0]}")
        for record_id, _, digest in records:
            calls[record_id] = call_from_row(by_hash[digest])
    else:
        import talea
        from talea.inference import TaleaPredictor
        from talea.io import SequenceRecord

        version = talea.__version__
        sequence_records = [SequenceRecord(record_id=record_id, sequence=sequence, sequence_sha256=digest)
                            for record_id, sequence, digest in records]
        with TaleaPredictor(device="mps", trust_remote_code=True, local_files_only=True) as predictor, \
                open(args.rows, "w") as handle:
            for row in predictor.iter_predictions(sequence_records):
                handle.write(json.dumps(row, sort_keys=True) + "\n")
                calls[row["record_id"]] = call_from_row(row)
    for record_id, numbers in residue_numbers.items():
        calls[record_id]["resi"] = numbers

    kept = {item["id"] for set_id, items in manifest["sets"].items() if set_id not in model_dirs for item in items}
    kept &= set(previous["calls"])
    if kept and previous.get("talea_version") != version:
        raise SystemExit(f"{calls_path} holds Talea {previous.get('talea_version')} calls; these are Talea {version}")
    merged = {**{key: previous["calls"][key] for key in previous["calls"] if key in kept}, **calls}
    payload = {"talea_version": version, "positions": "1-based inclusive", "calls": merged}
    calls_path.write_text(json.dumps(payload, separators=(",", ":")))
    with_regions = sum(1 for call in calls.values() if call["regions"])
    typed = sum(1 for call in calls.values() if "protein_type" in call)
    print(f"talea {version}: {len(calls)} proteins, {with_regions} with regions, "
          f"{sum(len(c['regions']) for c in calls.values())} regions, {typed} with types; "
          f"{len(kept)} calls kept from other sets")


if __name__ == "__main__":
    main()
