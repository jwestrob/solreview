"""Bench data for the Giants tab: one file per giant with its models, Talea's calls and tracks, plus its attention strip.

Inputs: the Talea run (talea_giants.py: rows and window captures), the structure bundles (encode_structures.py)
and the cohort table (TSV: record_id, length and an optional display name). For each giant writes into
BENCH/data/giants/:

- NN.png: uint32 header length | header JSON | model records | tracks, carried as the pixels of an 8-bit grayscale
  PNG (the artifact host serves images, not raw binary; the PNG's zlib stream compresses the bytes and its Adler-32
  checks them; the page parses the PNG itself and never decodes it as a picture). The header holds the id, display
  name, sequence, the models ([scheme, start, end, record offset, residues, mean pLDDT, pTM, ranking score, fold
  mode, substitutions]; positions 0-based half-open in the cohort sequence), Talea's windows and repeat regions
  (same frame) and where the tracks sit: the U-Net residue-state probability per residue and the helix and strand
  probes per 2 residues, all as bytes (value x 255).
- NN.jpg: the attention strip (build_giant_tracks.strip_image), layer 13 above layer 35.

and index.json listing the giants with their headline numbers.

usage: python3 package_giants.py RUN_DIR BUNDLE_DIR COHORT_TSV BENCH_DIR
"""
import csv
import json
import struct
import sys
import zlib
from pathlib import Path

import numpy as np

import build_giant_tracks as tracks

X_BIN = 8
JPEG_QUALITY = 80


def png_container(payload, width=4096):
    """An 8-bit grayscale PNG whose pixels are uint32 length (little-endian) + payload, zero-padded to whole rows."""
    data = struct.pack("<I", len(payload)) + payload
    rows = -(-len(data) // width)
    data += bytes(rows * width - len(data))
    scanlines = b"".join(b"\0" + data[row * width:(row + 1) * width] for row in range(rows))

    def chunk(kind, body):
        return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body) & 0xFFFFFFFF)

    header = struct.pack(">IIBBBBB", width, rows, 8, 0, 0, 0, 0)
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", header) + chunk(b"IDAT", zlib.compress(scanlines, 9))
            + chunk(b"IEND", b""))


def talea_summary(row):
    windows = [
        [w["start"], w["end"], round(w["talea_discovery_score"], 3), w["period"], w["wedge_start"], w["wedge_end"],
         w["registered_start"], w["registered_end"], w["n_units"], round(w["teacher_state_z"], 2)]
        for w in row["windowing"]["windows"]
    ]
    regions = []
    for region in row.get("repeat_regions") or []:
        units = [value for unit in region.get("units") or [] for value in (unit["start"], unit["end"])]
        regions.append([
            region["start"], region["end"], region["period"], 0 if region["source"] == "primary" else 1,
            region.get("window_index"), round(region.get("mean_state_probability") or 0.0, 3),
            round(region["wedge_score"], 2) if region.get("wedge_score") is not None else None, units,
        ])
    return {
        "best_window": row["windowing"]["best_window_index"],
        "score": round(row["talea_discovery_score"], 3),
        "windows": windows,
        "regions": regions,
    }


def main():
    run, bundle_dir, cohort_tsv, bench = (Path(value) for value in sys.argv[1:5])
    out = bench / "data" / "giants"
    out.mkdir(parents=True, exist_ok=True)
    rows = {json.loads(line)["record_id"]: json.loads(line) for line in open(run / "talea_giants_rows.jsonl")}
    bundles = json.loads((bundle_dir / "bundles.json").read_text())
    cohort = list(csv.DictReader(open(cohort_tsv), delimiter="\t"))
    index = []
    for entry in cohort:
        protein = entry["record_id"]
        if protein not in rows:
            print(f"skip {protein}: no Talea row yet", flush=True)
            continue
        row, bundle = rows[protein], bundles[protein]
        raw = zlib.decompress((bundle_dir / "bundles" / bundle["file"]).read_bytes())
        header_length = struct.unpack_from("<I", raw)[0]
        structure = json.loads(raw[4:4 + header_length])
        records = raw[4 + header_length:]
        length = len(structure["sequence"])
        assert length == row["length"] == int(entry["length"]), protein

        ends = {window["start"]: window["end"] for window in row["windowing"]["windows"]}
        starts = list(ends)
        cache = {}

        def capture(start):
            if start not in cache:
                cache.clear()
                cache[start] = np.load(run / "windows" / f"{protein}__{start + 1}-{ends[start]}.npz")
            return cache[start]

        strips = []
        for layer in (13, 35):
            band = tracks.blend(length, starts, lambda s: capture(s)[f"band{layer}"].astype(np.float32),
                                tracks.BAND_BIN)
            strips.append(tracks.strip_image(band, X_BIN))
        columns = min(strip.shape[1] for strip in strips)
        image = np.concatenate([strip[:, :columns] for strip in strips], axis=0)
        from PIL import Image
        stem = bundle["file"].removesuffix(".bin")
        Image.fromarray(tracks.magma(image)).save(out / f"{stem}.jpg", quality=JPEG_QUALITY)

        state = tracks.blend(length, starts, lambda s: capture(s)["state"].astype(np.float32))
        ss = tracks.blend(length, starts, lambda s: capture(s)["ss"].astype(np.float32).T)
        pairs = length // tracks.SS_BIN * tracks.SS_BIN
        ss = np.nanmean(ss[:2, :pairs].reshape(2, -1, tracks.SS_BIN), axis=2)
        as_bytes = lambda values: np.clip(np.rint(np.nan_to_num(values) * 255), 0, 255).astype(np.uint8).tobytes()
        state_bytes, ss_bytes = as_bytes(state), as_bytes(ss)

        header = {
            "id": protein,
            "name": entry.get("name") or protein,
            "length": length,
            "sequence": structure["sequence"],
            "models": structure["windows"],
            "talea": talea_summary(row),
            "tracks": {
                "state": [len(records), length],
                "ss": [len(records) + len(state_bytes), len(ss_bytes) // 2, tracks.SS_BIN],
            },
            "strip": {"columns": columns, "x_bin": X_BIN, "lag_step": tracks.LAG_STEP, "lags": tracks.LAGS,
                      "rows_per_layer": tracks.LAGS // tracks.LAG_STEP, "layers": [13, 35]},
        }
        encoded = json.dumps(header, separators=(",", ":")).encode()
        payload = struct.pack("<I", len(encoded)) + encoded + records + state_bytes + ss_bytes
        (out / f"{stem}.png").write_bytes(png_container(payload))
        talea = header["talea"]
        index.append({
            "id": protein,
            "name": header["name"],
            "length": length,
            "bundle": f"data/giants/{stem}.png",
            "strip": f"data/giants/{stem}.jpg",
            "models": len(structure["windows"]),
            "score": talea["score"],
            "regions": len(talea["regions"]),
            "windows": len(talea["windows"]),
            "x_residues": structure["sequence"].count("X"),
        })
        print(f"{stem} {header['name'][:40]:40s} {length:>6} aa {len(structure['windows']):>5} models "
              f"{len(talea['regions']):>3} regions {(out / (stem + '.png')).stat().st_size / 1e6:5.2f} MB "
              f"strip {(out / (stem + '.jpg')).stat().st_size / 1e6:4.2f} MB", flush=True)
    # talea_giants.py logs "talea <version>: ..." on its first line.
    version = (run / "talea_giants.log").read_text().split(":", 1)[0].removeprefix("talea ").strip()
    meta = {"talea_version": version, "proteins": index}
    (out / "index.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
