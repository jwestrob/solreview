"""Whole-protein tracks for the giants from Talea's window captures (talea_giants.py).

Talea scores a giant as 2,000-aa windows every 1,000 residues. Each window's capture holds its attention band
(residue pairs up to 300 apart, pair midpoints binned every 4 residues), its U-Net residue-state track and its
helix/strand/loop probes. This script blends overlapping windows with triangular weights that peak at each window's
centre (they sum to one where two windows overlap), so every residue takes most of its value from the window in
which it sits most centrally, and writes per protein:

- strip JPEG: the attention band of layers 13 and 35 stacked (each LAG_ROWS rows: lag 2 at the top to 300 at the
  bottom, in 2-residue steps; one column per X_BIN residues of pair midpoint), each lag row scaled from its median
  over the protein to its 99.7th percentile (strip_image) and coloured with magma;
- tracks: uint8 residue-state probability per residue and helix/strand probe probabilities per SS_BIN residues.

usage: python3 build_giant_tracks.py RUN_DIR OUT_DIR [--quality Q] [--xbin N]
"""
import argparse
import json
from pathlib import Path

import numpy as np

LAGS = 300
BAND_BIN = 4
LAG_STEP = 2
SS_BIN = 2


def window_weights(length):
    position = np.arange(length) + 0.5
    return np.minimum(position, length - position)


def blend(length, starts, values_of, bin_width=1):
    """Weighted mean over windows of per-position values (last axis: positions in bins of bin_width residues)."""
    columns = -(-length // bin_width)
    total = weight = None
    for start in starts:
        values = np.asarray(values_of(start), dtype=np.float64)
        offset = int(round(start / bin_width))
        n = min(values.shape[-1], columns - offset)
        taper = window_weights(values.shape[-1])[:n]
        part = values[..., :n]
        valid = np.isfinite(part)
        if total is None:
            total = np.zeros(values.shape[:-1] + (columns,))
            weight = np.zeros_like(total)
        total[..., offset:offset + n] += np.where(valid, part, 0.0) * taper
        weight[..., offset:offset + n] += valid * taper
    with np.errstate(invalid="ignore", divide="ignore"):
        return total / weight


def magma(values):
    from matplotlib import colormaps

    return (colormaps["magma"](np.clip(values, 0.0, 1.0))[..., :3] * 255).astype(np.uint8)


def strip_image(band, x_bin):
    """Lag rows 1..LAGS (4-residue columns) -> LAGS/LAG_STEP rows and x_bin-residue columns, scaled 0..1.

    Attention falls off with sequence separation, so each lag row is scaled on its own: 0 at the row's median over
    the protein, 1 at its 99.7th percentile. Periodic stripes then stay visible out to lag 300.
    """
    lag_rows = band.reshape(LAGS // LAG_STEP, LAG_STEP, band.shape[1])
    merge = x_bin // BAND_BIN
    columns = band.shape[1] // merge * merge
    merged = np.nanmean(lag_rows[:, :, :columns].reshape(LAGS // LAG_STEP, LAG_STEP, -1, merge), axis=(1, 3))
    median = np.nanmedian(merged, axis=1, keepdims=True)
    high = np.nanpercentile(merged, 99.7, axis=1, keepdims=True)
    scaled = (merged - median) / np.maximum(high - median, 1e-12)
    return np.nan_to_num(np.clip(scaled, 0.0, 1.0), nan=0.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("out")
    parser.add_argument("--quality", type=int, default=80)
    parser.add_argument("--xbin", type=int, default=8)
    parser.add_argument("--only", nargs="*")
    args = parser.parse_args()
    from PIL import Image

    run, out = Path(args.run), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    summary = {}
    for line in open(run / "talea_giants_rows.jsonl"):
        row = json.loads(line)
        protein, length = row["record_id"], row["length"]
        if args.only and protein not in args.only:
            continue
        starts = [window["start"] for window in row["windowing"]["windows"]]
        cache = {}

        def capture(start):
            if start not in cache:
                cache.clear()
                cache[start] = np.load(run / "windows" / f"{protein}__{start + 1}-{start + 2000}.npz")
            return cache[start]

        layers = []
        for layer in (13, 35):
            band = blend(length, starts, lambda s: capture(s)[f"band{layer}"].astype(np.float32), BAND_BIN)
            layers.append(strip_image(band, args.xbin))
        state = blend(length, starts, lambda s: capture(s)["state"].astype(np.float32))
        ss = blend(length, starts, lambda s: capture(s)["ss"].astype(np.float32).T)
        columns = min(layer.shape[1] for layer in layers)
        rows = np.concatenate([layer[:, :columns] for layer in layers], axis=0)
        Image.fromarray(magma(rows)).save(out / f"{protein}.jpg", quality=args.quality)
        ss_binned = np.nanmean(ss[:2, : length // SS_BIN * SS_BIN].reshape(2, -1, SS_BIN), axis=2)
        np.savez(out / f"{protein}_tracks.npz",
                 state=np.clip(np.rint(np.nan_to_num(state) * 255), 0, 255).astype(np.uint8),
                 ss=np.clip(np.rint(np.nan_to_num(ss_binned) * 255), 0, 255).astype(np.uint8))
        size = (out / f"{protein}.jpg").stat().st_size
        summary[protein] = {"length": length, "strip_columns": columns, "strip_bytes": size}
        print(protein[:50], length, columns, f"{size / 1e3:.0f} kB", flush=True)
    (out / "tracks_summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
