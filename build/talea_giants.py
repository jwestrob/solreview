"""Talea on long proteins in 2,000-aa windows, keeping each window's attention for the review bench.

Runs the installed Talea (MPS) in window mode over every protein of FASTA and writes one row per protein to
OUT_DIR/talea_giants_rows.jsonl. Each window's postprocessing task also saves, in the worker, OUT_DIR/windows/
<protein>__<start>-<end>.npz with float16 arrays:

- band<layer>: the APC attention band, residue pairs up to LAGS apart (row d-1 holds pairs d apart, the mean of
  both triangles), columns binning pair midpoints every BIN residues of the window;
- square<layer>: the whole window matrix as SQUARE x SQUARE block means;
- state: the U-Net residue-state track; ss: the layer-35 helix/strand/loop probes (residues x 3).

Talea's contract is the 20 canonical residues; here X (at assembly gaps, for example) goes to ESM++ as its unknown
token, and only the windows that contain it are affected.

usage: python3 talea_giants.py FASTA OUT_DIR
"""
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

LAGS = 300
BIN = 4
SQUARE = 400
WINDOW_DIR = "TALEA_GIANTS_WINDOW_DIR"


def attention_band(matrix, lags, bin_width):
    length = matrix.shape[0]
    columns = -(-length // bin_width)
    band = np.full((lags, columns), np.nan, dtype=np.float32)
    for lag in range(1, min(lags, length - 1) + 1):
        values = 0.5 * (np.diagonal(matrix, lag) + np.diagonal(matrix, -lag))
        # Twice the pair midpoint i + lag/2, so odd lags bin without rounding.
        column = (2 * np.arange(length - lag) + lag) // (2 * bin_width)
        sums = np.bincount(column, values, minlength=columns)
        counts = np.bincount(column, minlength=columns)
        band[lag - 1] = np.where(counts > 0, sums / np.maximum(counts, 1), np.nan)
    return band


def block_mean(matrix, size):
    edges = np.linspace(0, matrix.shape[0], size + 1).astype(int)
    starts, widths = edges[:-1], np.diff(edges)
    rows = np.add.reduceat(matrix, starts, axis=0) / widths[:, None]
    return np.add.reduceat(rows, starts, axis=1) / widths[None, :]


def capture_run_record(payload):
    """Save the window's attention band, square, state and probes, then postprocess it as Talea would."""
    from talea.postprocess_pool import run_record

    metadata, matrices, _, probabilities, state = payload
    arrays = {}
    for layer, matrix in matrices.items():
        values = np.asarray(matrix, dtype=np.float32)
        arrays[f"band{layer}"] = attention_band(values, LAGS, BIN).astype(np.float16)
        arrays[f"square{layer}"] = block_mean(values, min(SQUARE, values.shape[0])).astype(np.float16)
    arrays["state"] = np.asarray(state, dtype=np.float16)
    if probabilities is not None:
        arrays["ss"] = np.asarray(probabilities, dtype=np.float16)
    name = metadata["record_id"].replace("/", "__")
    path = Path(os.environ[WINDOW_DIR]) / f"{name}.npz"
    np.savez(path.with_suffix(".tmp.npz"), **arrays)
    path.with_suffix(".tmp.npz").replace(path)
    return run_record(payload)


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


def main():
    fasta, out = Path(sys.argv[1]), Path(sys.argv[2])
    (out / "windows").mkdir(parents=True, exist_ok=True)
    os.environ[WINDOW_DIR] = str(out / "windows")

    import talea
    import talea.inference as inference
    from talea.io import SequenceRecord

    class ExplorationPredictor(inference.TaleaPredictor):
        def _validate_records(self, records):
            allowed = frozenset("ACDEFGHIKLMNPQRSTVWYX")
            for record in records:
                unexpected = sorted(set(record.sequence) - allowed)
                if unexpected or record.length > self.maximum_length:
                    raise ValueError(f"{record.record_id}: {''.join(unexpected)} / {record.length} aa")
                if hashlib.sha256(record.sequence.encode()).hexdigest() != record.sequence_sha256:
                    raise ValueError(f"sequence checksum mismatch for {record.record_id}")

    inference.run_record = capture_run_record
    records = [
        SequenceRecord(record_id=name, sequence=sequence,
                       sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest())
        for name, sequence in read_fasta(fasta).items()
    ]
    started = time.time()
    print(f"talea {talea.__version__}: {len(records)} proteins, {sum(r.length for r in records):,} aa", flush=True)
    with ExplorationPredictor(device="mps", trust_remote_code=True, local_files_only=True) as predictor, \
            open(out / "talea_giants_rows.jsonl", "w") as handle:
        for done, row in enumerate(predictor.iter_windowed_predictions(records), 1):
            handle.write(json.dumps(row, sort_keys=True) + "\n")
            handle.flush()
            print(f"{done}/{len(records)} {row['record_id']} {row['length']} aa, "
                  f"{row['windowing']['window_count']} windows, best {row['talea_discovery_score']:.2f}, "
                  f"{len(row.get('repeat_regions') or [])} regions, {time.time() - started:.0f} s", flush=True)


if __name__ == "__main__":
    main()
