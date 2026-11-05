# Building a dataset

`BENCH` below is the folder that holds `data/`; `BENCH=.` writes into this repo's `data/`. The scripts need numpy,
Pillow and matplotlib. The steps that run Talea also need Talea and its weights, either installed or on `PYTHONPATH`.
They run on Apple-silicon MPS.

## Review sets

Put the PDB models for each set in their own folder. Blinded file names such as `review_01.pdb` keep the review
blind.

```sh
PYTHONPATH=<Talea tree> python3 build/prepare_review_data.py BENCH census=models/census bfvd=models/bfvd
PYTHONPATH=<Talea tree> python3 build/talea_calls.py BENCH rows.jsonl census=models/census bfvd=models/bfvd
```

- `prepare_review_data.py` writes the gzipped models and attention maps to `BENCH/data/<set>/`, plus
  `BENCH/data/manifest.json`. Talea's scores stay out of the page. A protein longer than 2,046 aa is scored in
  Talea's 2,000-aa windows, and its map shows the window Talea reports.
- `talea_calls.py` writes Talea's regions to `BENCH/data/talea_calls.json` and its full rows to `rows.jsonl`. Rows
  that carry Talea's type head add each region's type and the whole protein's.
- Both scripts rebuild the sets named on the command line and keep the others.

To reuse rows from an earlier Talea run, such as a windowed run with region types, pass them to `talea_calls.py`
with `--from-rows`. The rows match the models by sequence, so Talea does not run again:

```sh
python3 build/talea_calls.py --from-rows --talea-version 0.3.2 BENCH talea_rows.jsonl dpann=models/dpann
```

Each call records the Talea version that wrote it, so rebuilt sets can run a newer Talea than the sets kept.

## Giants

1. **Run Talea** in window mode, saving each window's attention band, residue state and secondary-structure probes.
   The script wraps Talea's postprocessing worker and was written against Talea 0.3.1. `package_giants.py` reads
   Talea's version from the first line of the log.

   ```sh
   mkdir -p run && python3 build/talea_giants.py giants.faa run > run/talea_giants.log
   ```

2. **Pack the folded models.** List them in a TSV with the columns `protein`, `scheme`, `start` and `path`. The
   columns `mode`, `ptm`, `ranking` and `offsets` are optional (see the script's docstring). Schemes are named
   `w<size>_s<step>`. The `w2000_s1000` models line up with Talea's windows and feed the ranked list.

   ```sh
   python3 build/encode_structures.py giants.faa models.tsv bundles
   ```

3. **Package each protein** into its bundle and attention strip. `giants.tsv` has the columns `record_id`, `length`
   and an optional display `name`.

   ```sh
   python3 build/package_giants.py run bundles giants.tsv BENCH
   ```

4. **Add optional overlays** from your own annotation: `BENCH/data/giants/hits.json` and `window_gate.json`. Their
   formats are in [docs/data-format.md](../docs/data-format.md).

5. **Rank the 2,000-aa models** by Talea window score. Run this after the overlays: it takes repeat regions from
   `window_gate.json` and hits from `hits.json` when they exist. `--exclude` takes a TSV of spans to leave out
   (`protein`, `start`, `end`, 0-based), such as models already labelled.

   ```sh
   python3 build/build_ranked_review.py run bundles BENCH/data --exclude labelled.tsv
   ```

`backbone_codec.py` holds the 5-byte-per-residue model encoding, and `build_giant_tracks.py` holds the window
blending and the strip image. The page decodes both with the same arithmetic.
