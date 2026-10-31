"""Data for the browser structure review: gzipped models, attention maps and a manifest.

For each blinded review set, reads the sequence from every model's C-alpha records, runs Talea on MPS in thread
mode and renders the layer-13 and layer-35 APC attention matrices Talea works from (side by side, magma, clipped
at the 1st and 99.5th percentiles), then gzips the model. Talea's scores are discarded: the review stays blind.

usage: PYTHONPATH=<Talea tree> python3 prepare_review_data.py OUT_DIR SET_ID=MODEL_DIR [SET_ID=MODEL_DIR ...]
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
PANEL = 480
GAP = 8


def model_sequence(path):
    residues = []
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            residues.append(THREE_TO_ONE[line[17:20].strip()])
    return "".join(residues)


def render(matrices, layers, out_path):
    from matplotlib import colormaps
    from PIL import Image

    cmap = colormaps["magma"]
    panels = []
    for layer in layers:
        values = np.asarray(matrices[layer], dtype=np.float32)
        length = values.shape[0]
        off = values[~np.eye(length, dtype=bool)]
        low, high = np.percentile(off, [1.0, 99.5])
        scaled = np.clip((values - low) / max(high - low, 1e-12), 0.0, 1.0)
        size = min(length, PANEL)
        image = Image.fromarray(scaled.astype(np.float32), mode="F").resize((size, size), Image.BOX)
        panels.append((cmap(np.asarray(image))[..., :3] * 255).astype(np.uint8))
    size = panels[0].shape[0]
    canvas = np.empty((size, 2 * size + GAP, 3), dtype=np.uint8)
    canvas[:] = (np.asarray(cmap(0.0)[:3]) * 255).astype(np.uint8)
    canvas[:, :size] = panels[0]
    canvas[:, size + GAP :] = panels[1]
    Image.fromarray(canvas).save(out_path, quality=88)


def main():
    import talea.inference as inference
    from talea.io import SequenceRecord

    out = Path(sys.argv[1])
    sets = [argument.split("=", 1) for argument in sys.argv[2:]]
    records, items = [], {}
    for set_id, model_dir in sets:
        (out / "data" / set_id).mkdir(parents=True, exist_ok=True)
        for path in sorted(Path(model_dir).glob("*.pdb")):
            sequence = model_sequence(path)
            record_id = f"{set_id}-{path.stem}"
            records.append(SequenceRecord(record_id=record_id, sequence=sequence,
                                          sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest()))
            items[record_id] = {
                "id": record_id, "set": set_id, "file": path.name, "length": len(sequence),
                "pdb": f"data/{set_id}/{path.stem}.pdb.gz", "attention": f"data/{set_id}/{path.stem}_attention.jpg",
            }
            with open(path, "rb") as source, gzip.open(out / items[record_id]["pdb"], "wb", compresslevel=9) as target:
                target.write(source.read())

    original = inference.process_attention_record
    layers = None

    def capture(task):
        record_id = task["metadata"]["record_id"]
        render(task["matrices"], layers, out / items[record_id]["attention"])
        return original(task)

    inference.process_attention_record = capture
    predictor = inference.TaleaPredictor(device="mps", trust_remote_code=True, local_files_only=True)
    layers = tuple(predictor.expected_layers)
    done = 0
    for _ in predictor.iter_predictions(records, postprocess_mode="thread"):
        done += 1
    manifest = {"layers": list(layers), "sets": {}}
    for set_id, _ in sets:
        manifest["sets"][set_id] = [items[key] for key in sorted(items) if items[key]["set"] == set_id]
    (out / "data" / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"{done} proteins; layers {layers}; " + ", ".join(f"{k}: {len(v)}" for k, v in manifest["sets"].items()))


if __name__ == "__main__":
    main()
