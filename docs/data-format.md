# Review data

The page reads everything from `data/` (or the folder passed to `serve.py --data`). Paths inside the JSON files are
relative to the page, so they start with `data/`. A dataset holds review sets, giants, or both. Unless noted,
positions are 0-based and half-open (`start` inclusive, `end` exclusive), as in Talea's rows.

```
data/
  sets.json              optional: tab titles and descriptions
  manifest.json          review sets: one model per protein
  talea_calls.json       Talea's repeat regions for the review sets
  review_hits.json       optional: annotation hits for review-set models (build/review_hits.py)
  <set>/                 models and attention maps per review set
  giants/
    index.json           the long proteins
    NN.png, NN.jpg       per protein: bundle (models, Talea calls, tracks) and attention strip
    ranked.json          optional: 2,000-aa models by Talea window score
    hits.json            optional: Pfam, CDD and Foldseek hits
    window_gate.json     optional: an alternative set of repeat regions to compare with Talea's
  labels/
    <reviewer>.json      written by serve.py: one reviewer's labels, {label id: label} (README.md)
```

## sets.json

```json
{"census": {"title": "Census hits", "about": "One line shown under the tabs."}, "giants": {"title": "Giants", "about": "..."}}
```

Tabs follow this file's order, then any other review set in `manifest.json`, then Giants. A set without an entry
shows its id as the title. An optional `"order": "gap"` opens a review set ordered by Talea's whole-protein model z
minus its residue model z, highest first (needs `scores` in `talea_calls.json`); the default is the manifest order.
`"manifest_order": ["By family size", "Largest family first"]` renames the manifest-order button and its caption
(default `["By score", "Highest Talea score first"]`).

## Review sets

`manifest.json` (written by `build/prepare_review_data.py`):

```json
{"layers": [13, 35],
 "sets": {"census": [{"id": "census-review_01", "set": "census", "file": "review_01.pdb", "length": 676,
                      "pdb": "data/census/review_01.pdb.gz", "attention": "data/census/review_01_attention.jpg"}]}}
```

- `id` is the label id, `<set>-<file stem>`. Name files by the protein's own accession (`MCX8197437.1.pdb`), not by
  rank or set position, so the list shows an identifier people can look up.
- `pdb` is a PDB file, plain or gzipped, with pLDDT in the B-factor column.
- `attention` is a JPEG of Talea's two attention maps side by side: layer 13 on the left and layer 35 on the right.
  Each panel is min(L, 480) px square with an 8 px gap between them, coloured with magma and clipped at the 1st and
  99.5th percentiles of the off-diagonal values.
- `attention_start` and `attention_end` (optional) mark a map that covers part of the protein. Talea reads at most
  2,046 residues, so a longer protein is scored in 2,000-aa windows every 1,000 residues; its map shows the window
  Talea reports, and L in the panel size is that window's length.
- `organism`, `lineage` and `genome` (optional) describe where the protein comes from. The list shows `organism`
  under the name (keep it short: phylum and order, or the species and locus tag). The viewer header shows `lineage`
  (e.g. the GTDB ranks) and `genome` (the assembly or bin accession).
- `alias` (optional) is an earlier name for the item; the page shows it on hover as "formerly …".

`talea_calls.json` (written by `build/talea_calls.py`):

```json
{"talea_versions": ["0.3.1", "0.3.2"], "positions": "1-based inclusive",
 "calls": {"census-review_06": {"talea_version": "0.3.1",
                                "regions": [{"start": 135, "end": 217, "period": 32.0, "source": "primary",
                                             "units": [{"start": 135, "end": 148}]}],
                                "resi": [5, 6, 7]},
           "dpann-dpann_001": {"talea_version": "0.3.2",
                               "regions": [{"start": 1, "end": 136, "period": 10.0, "source": "primary",
                                            "units": [{"start": 1, "end": 4}],
                                            "type": {"alpha": 0.009, "alpha_beta": 0.156, "beta": 0.729,
                                                     "beta_propeller_or_beta_sandwich": 0.005, "non_solenoid": 0.101}},
                                           {"start": 160, "end": 214, "period": null, "source": "state", "units": []}],
                               "protein_type": {"alpha": 0.009, "alpha_beta": 0.156, "beta": 0.729,
                                                "beta_propeller_or_beta_sandwich": 0.005, "non_solenoid": 0.101}}}}
```

- Positions here are 1-based and inclusive, counted along the modelled sequence.
- `talea_version` is the Talea that wrote the call, and `talea_versions` lists every version in the file. The page
  names each protein's version above its regions. An older file has one top-level `talea_version` for every call.
- `source` is Talea's: `primary`, `secondary` for the further wedge regions, or `state` (Talea 0.3.2 and later) for
  a region from the residue-state track alone. A state region marks a stretch the residue model flags where the
  attention search found no period, so it has `period: null` and no units. The page shows `secondary` regions as
  further regions. It draws a state region as a dotted outline on the maps and one flat colour on the structure, and
  lists it as "residue model only · no period found".
- `resi` appears only for a model numbered other than 1..L, and gives the residue number at each position.
- `scores` (optional): Talea's score and its two halves, `{"talea": 2.41, "residue_z": 0.87, "whole_protein_z": 3.92}`
  (`talea_discovery_score`, `teacher_state_z`, `global_z`). Sets with scores show them beside Talea's regions and can be
  listed by `whole_protein_z - residue_z`.
- `type` and `protein_type` (optional) are Talea's five-class type head: the probabilities of alpha solenoid, beta
  solenoid, alpha/beta solenoid, beta propeller or beta sandwich, and not a solenoid. `type` averages the head over
  the region and `protein_type` over the whole protein. They come from Talea rows that carry
  `region_type_probability` on each region and `predicted_teacher_mean_class_probability`. The page shows the
  primary region's type with Talea's regions, as an unverified prediction.
- For a protein Talea scores in windows, `protein_type` covers the window Talea reports. Each region keeps the type
  from the window that called it, averaged over that window's span of the region, including where Talea extends the
  region with calls from other windows.

## Giants

`giants/index.json` (written by `build/package_giants.py`):

```json
{"talea_version": "0.3.1",
 "proteins": [{"id": "<record id>", "name": "<display name>", "length": 87279, "bundle": "data/giants/00.png",
               "strip": "data/giants/00.jpg", "models": 783, "score": 2.836, "regions": 46, "windows": 87,
               "x_residues": 0}]}
```

The page lists the proteins longest first. The bundle's file stem (`00`) is the protein's key in `ranked.json` and
in giant label ids.

### NN.png: the bundle

The bundle is a binary payload carried as the pixels of an 8-bit grayscale PNG. The page parses the PNG byte by byte
and never decodes it as a picture, so any web host that serves images can serve it. The pixel bytes are a uint32
little-endian length `N`, then `N` payload bytes, then zero padding to the end of the row. The payload is:

```
uint32 LE header length H | header JSON (H bytes, UTF-8) | body
```

The header holds:

- `id`, `name`, `length`, `sequence`.
- `models`: one row per folded model, `[scheme, start, end, body offset, residues, mean pLDDT, pTM, ranking score,
  mode, substitutions]`.
  - `scheme` is `w<size>_s<step>`, for example `w2000_s1000`. The page offers one button per window size.
  - `body offset` is where the model's record starts in the body.
  - `mode` is free text shown under the model; `substitutions` counts residues that differ from the cohort sequence.
- `talea`, with these fields:
  - `best_window`: the index of the best window.
  - `score`: Talea's discovery score for the protein.
  - `windows`: one row per window, `[start, end, score, period, wedge start, wedge end, registered start, registered
    end, units, teacher state z]`.
  - `regions`: one row per region, `[start, end, period, further (0 primary, 1 further), window index, mean state
    probability, wedge score, [unit start, unit end, ...]]`.
- `tracks`:
  - `state: [offset, n]` is the U-Net residue-state probability per residue, one byte each (value × 255).
  - `ss: [offset, n, bin]` is the helix probe, then the strand probe, each `n` bytes averaging `bin` residues.
- `strip`: the layout of the attention strip, `{columns, x_bin, lag_step, lags, rows_per_layer, layers}`.

The body holds the model records back to back, then the state bytes, then the probe bytes. A model record is
`build/backbone_codec.py`'s compact backbone, 5 bytes per residue:

- uint32 `n`, the number of residues;
- int16[3], the first C-alpha on a 0.1 Å grid;
- uint16 `k`, then `k` × (uint32 step index, int16[3] step) for C-alpha steps longer than 12.7 Å;
- int8[3(n−1)] C-alpha steps, with 0 at the steps listed above;
- uint8[n−1] peptide-plane angles;
- uint8[n] pLDDT.

The page rebuilds N, C and O from the trace and the plane angles.

### NN.jpg: the attention strip

The strip shows Talea's attention band along the whole protein: layer 13 above, layer 35 below. Each layer takes
`rows_per_layer` rows, from residue pairs `lag_step` apart at the top to `lags` apart at the bottom. Columns step
along the pair midpoint every `x_bin` residues. Overlapping windows blend with triangular weights centred on each
window. Attention falls off with separation, so each row is scaled on its own, from its median over the protein
to its 99.7th percentile, and coloured with magma.

### ranked.json (optional)

`ranked.json` drives the Giants list's "By score" view and is written by `build/build_ranked_review.py`:

```json
{"rule": "...", "items": [["00", 1000, 3000, 2.512, [[1210, 2870, 34.5, 0]], "CDD", "TPR", [["31", 5000, 7000, 2.31]]]]}
```

Each item is:

- the protein key;
- the 2,000-aa model's start and end, and its window score;
- the repeat regions inside the model, `[start, end, period, new]`;
- the hit source (`Pfam`, `CDD`, `Pfam/CDD`, `Foldseek` or empty) and the hit text;
- near-identical copies of the model in other proteins, `[key, start, end, score]`.

Items are sorted by score, highest first.

### review_hits.json (optional)

```json
{"sources": "one line on where the hits come from",
 "items": {"dpann_gap-gap_004": [["Pfam", 20, 87, "LpxD", "UDP-3-O-[3-hydroxymyristoyl] glucosamine N-acyltransferase", 9.9e-21],
                                 ["KOfam", null, null, "K02536", "lpxD", 1e-30]]}}
```

Each hit is `[source, start, end, family, description, E-value]`, with start and end 0-based end-exclusive along the
modelled sequence, or both null for a whole-protein call. The sidebar lists the hits for the model on screen, as
observed annotations. Written by `build/review_hits.py` from a TSV.

### hits.json (optional)

```json
{"sources": "one line on where the hits come from",
 "giants": {"<record id>": {"domains": [["Pfam", 120, 410, "TPR_1", "Tetratricopeptide repeat", 1e-8, 45.2]],
                            "foldseek": [[2000, 2410, "AF-Q9XXX-F1", 0, 1e-12, 310, 0.31, 3]]}},
 "descriptions": ["description of target 0"]}
```

- `domains`: `[source, start, end, family, description, E-value, score]`.
- `foldseek`: `[start, end, target, description index, E-value, bits, identity, window]`. An AlphaFold DB target
  becomes a link.

The sidebar lists the hits that overlap the model on screen.

### window_gate.json (optional)

The overlay is an alternative set of repeat regions, such as a different gate on Talea's window merging. With it,
the track toggles between Talea's regions and the overlay's:

```json
{"rule": "one line describing the alternative",
 "giants": {"<record id>": {"regions": [[374, 488, 28.5, 1, 0, 0.144, 18.95, [374, 401, 401, 488], 1]],
                            "candidates": [[374, 488, 28.5, 0, 0, 0.144, 18.95, 2.109, 0, 1]]}}}
```

- `regions` use the bundle's region layout, followed by `new` (1 when the region is mostly outside Talea's own
  regions).
- `candidates` are every region either set could hold, `[start, end, period, kind (0 window primary, 1 peeled call),
  window, mean state, wedge score, window score, kept by Talea, kept by the overlay]`. The track outlines the ones
  the shown set drops.
