"""Talea's region calls for the review bench: every repeat region with its units, positions 1-based inclusive.

Runs the Talea tree on PYTHONPATH (MPS) over every model in the bench manifest, writes the full rows to ROWS_JSONL
(kept off the page) and, into the bench data, talea_calls.json with only where Talea places repeats: region
extents, period, primary or further, and units, as sequence positions; a model numbered other than 1..L also gets
`resi`, the residue number at each position. Scores stay out of the page.

usage: PYTHONPATH=<Talea tree> python3 talea_calls.py BENCH_DIR ROWS_JSONL SET_ID=MODEL_DIR [SET_ID=MODEL_DIR ...]
"""
import hashlib
import json
import sys
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


def main():
    import talea
    from talea.inference import TaleaPredictor
    from talea.io import SequenceRecord

    bench, rows_path = Path(sys.argv[1]), Path(sys.argv[2])
    model_dirs = dict(argument.split("=", 1) for argument in sys.argv[3:])
    manifest = json.loads((bench / "data" / "manifest.json").read_text())
    records, residue_numbers = [], {}
    for set_id, items in manifest["sets"].items():
        for item in items:
            numbers, sequence = model_residues(Path(model_dirs[set_id]) / item["file"])
            # Talea positions index the modelled sequence; a model numbered other than 1..L also gets the
            # residue number of every position, so the page can show regions in the model's own numbering.
            if numbers != list(range(1, len(sequence) + 1)):
                residue_numbers[item["id"]] = numbers
            records.append(SequenceRecord(record_id=item["id"], sequence=sequence,
                                          sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest()))
    calls = {}
    with TaleaPredictor(device="mps", trust_remote_code=True, local_files_only=True) as predictor, \
            open(rows_path, "w") as handle:
        for row in predictor.iter_predictions(records):
            handle.write(json.dumps(row, sort_keys=True) + "\n")
            regions = []
            for region in row.get("repeat_regions") or []:
                regions.append({
                    **one_based(region),
                    "period": round(float(region["period"]), 1),
                    "source": region["source"],
                    "units": [one_based(unit) for unit in region.get("units") or []],
                })
            calls[row["record_id"]] = {"regions": regions}
            if row["record_id"] in residue_numbers:
                calls[row["record_id"]]["resi"] = residue_numbers[row["record_id"]]
    payload = {"talea_version": talea.__version__, "positions": "1-based inclusive", "calls": calls}
    (bench / "data" / "talea_calls.json").write_text(json.dumps(payload, separators=(",", ":")))
    with_regions = sum(1 for call in calls.values() if call["regions"])
    print(f"talea {talea.__version__}: {len(calls)} proteins, {with_regions} with regions, "
          f"{sum(len(c['regions']) for c in calls.values())} regions")


if __name__ == "__main__":
    main()
