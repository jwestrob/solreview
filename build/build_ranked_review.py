"""The Giants tab's ranked review list: 2,000-aa models in order of Talea window score, highest first.

Talea's window mode scores the same 2,000/1,000 windows the bench folds, so each window maps to the model holding
its centre. Runs of sequential windows collapse by greedy suppression along each protein: windows are taken in order
of score and a window enters only when it shares no residue with a window already taken. Near-identical copies then
fold together across proteins: a model sharing at least half of its 6-mers with a higher-ranked entry (about 90%
identity or more) becomes a copy of that entry. Spans listed in EXCLUDE_TSV (columns protein, start, end; 0-based,
half-open), such as models already labelled, leave the list together with their near-identical copies.

Each entry carries the repeat regions in its model (half of the region or half of the model inside): the overlay's
regions when DATA_DIR/giants/window_gate.json exists, else Talea's own. With DATA_DIR/giants/hits.json it also
carries the Pfam or CDD domains under those regions, else the strongest named Foldseek hit there (observed hits;
nothing is named here).

Writes DATA_DIR/giants/ranked.json: {"rule": ..., "items": [[key, model start0, model end, window score,
[[start0, end, period, new]], hit source, hit text, [[copy key, copy start0, copy end, copy score]]]]}, items
sorted by score, key being the protein's bundle name in DATA_DIR/giants/index.json.

usage: python3 build_ranked_review.py RUN_DIR BUNDLE_DIR DATA_DIR [--exclude EXCLUDE_TSV]
"""
import argparse
import csv
import json
import re
import struct
import zlib
from pathlib import Path

SCHEME = "w2000_s1000"
UNNAMED = re.compile(r"uncharacteri[sz]ed|hypothetical|unknown function|^$", re.I)
K = 6
SHARED_KMERS = 0.5


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def read_run(run_dir):
    """Talea's windows [(start0, end, score)] and repeat regions [[start0, end, period, 0]] per protein."""
    windows, regions = {}, {}
    for line in open(Path(run_dir) / "talea_giants_rows.jsonl"):
        row = json.loads(line)
        windows[row["record_id"]] = sorted((w["start"], w["end"], float(w["talea_discovery_score"]))
                                           for w in row["windowing"]["windows"])
        regions[row["record_id"]] = [[r["start"], r["end"], r["period"], 0] for r in row.get("repeat_regions") or []]
    return windows, regions


def read_bundles(bundle_dir):
    """Each protein's sequence and its SCHEME models [(start0, end)]."""
    manifest = json.load(open(Path(bundle_dir) / "bundles.json"))
    sequences, models = {}, {}
    for protein, entry in manifest.items():
        raw = zlib.decompress((Path(bundle_dir) / "bundles" / entry["file"]).read_bytes())
        header = json.loads(raw[4:4 + struct.unpack("<I", raw[:4])[0]])
        sequences[protein] = header["sequence"]
        models[protein] = sorted((w[1], w[2]) for w in header["windows"] if w[0] == SCHEME)
    return sequences, models


def model_at(models, position):
    """The model holding ``position`` nearest its centre (the page's modelIndexAt), or None."""
    inside = [m for m in models if m[0] <= position < m[1]]
    return min(inside, key=lambda m: abs((m[0] + m[1]) / 2 - position)) if inside else None


def hit_tag(calls, model, hits, descriptions):
    """(source, text) for the calls inside a model, or for the whole model when it holds none."""
    spans = [(c[0], c[1]) for c in calls] or [model]
    near = lambda a, b: any(overlap(a, b, s, e) >= 0.5 * min(b - a, e - s) for s, e in spans)
    domains = sorted((h for h in hits["domains"] if near(h[1], h[2])), key=lambda h: h[5] if h[5] is not None else 1.0)
    if domains:
        names = list(dict.fromkeys(h[3] for h in domains))[:2]
        sources = {h[0] for h in domains if h[3] in names}
        return ("Pfam" if sources == {"Pfam"} else "CDD" if sources == {"CDD"} else "Pfam/CDD"), " · ".join(names)
    fold = sorted((h for h in hits["foldseek"] if near(h[0], h[1]) and not UNNAMED.search(descriptions[h[3]])),
                  key=lambda h: -h[5])
    return ("Foldseek", descriptions[fold[0][3]]) if fold else ("", "")


def kmers(sequence):
    return {sequence[i:i + K] for i in range(len(sequence) - K + 1)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("run_dir", type=Path, help="talea_giants.py output (talea_giants_rows.jsonl)")
    parser.add_argument("bundle_dir", type=Path, help="encode_structures.py output (bundles.json, bundles/)")
    parser.add_argument("data_dir", type=Path, help="the bench data folder (giants/index.json)")
    parser.add_argument("--exclude", type=Path, help="TSV of spans to leave out: protein, start, end (0-based)")
    args = parser.parse_args()
    scores, talea_regions = read_run(args.run_dir)
    sequences, models = read_bundles(args.bundle_dir)
    giants = args.data_dir / "giants"
    index = json.load(open(giants / "index.json"))["proteins"]
    key = {p["id"]: Path(p["bundle"]).stem for p in index}
    overlay = json.load(open(giants / "window_gate.json"))["giants"] if (giants / "window_gate.json").exists() else None
    hits = json.load(open(giants / "hits.json")) if (giants / "hits.json").exists() else None
    excluded = {}
    for row in csv.DictReader(open(args.exclude), delimiter="\t") if args.exclude else []:
        excluded.setdefault(row["protein"], []).append((int(row["start"]), int(row["end"])))

    candidates, blockers, left_out = [], [], {"excluded": 0, "no model": 0, "copies of excluded": 0}
    for protein, windows in scores.items():
        kept = []
        for start, end, score in sorted(windows, key=lambda w: (-w[2], w[0])):
            if all(overlap(start, end, a, b) == 0 for a, b, _ in kept):
                kept.append((start, end, score))
        for start, end, score in kept:
            model = model_at(models.get(protein, []), (start + end) / 2)
            if model is None:
                left_out["no model"] += 1
            elif any(overlap(*model, a, b) > 0.5 * (model[1] - model[0]) for a, b in excluded.get(protein, [])):
                left_out["excluded"] += 1
                blockers.append(kmers(sequences[protein][model[0]:model[1]]))
            else:
                candidates.append((score, protein, model))
    for protein, spans in excluded.items():
        blockers.extend(kmers(sequences[protein][a:b]) for a, b in spans)

    def copies(words, pool):
        shared = lambda other: len(words & other) >= SHARED_KMERS * min(len(words), len(other))
        return next((i for i, other in enumerate(pool) if shared(other)), None)

    items, taken = [], []
    for score, protein, model in sorted(candidates, key=lambda c: (-c[0], key[c[1]], c[2][0])):
        words = kmers(sequences[protein][model[0]:model[1]])
        if copies(words, blockers) is not None:
            left_out["copies of excluded"] += 1
            continue
        twin = copies(words, taken)
        if twin is not None:
            items[twin][7].append([key[protein], model[0], model[1], round(score, 3)])
            continue
        regions = [[r[0], r[1], r[2], r[8]] for r in overlay[protein]["regions"]] if overlay else talea_regions[protein]
        calls = [r for r in regions if overlap(r[0], r[1], *model) >= 0.5 * min(r[1] - r[0], model[1] - model[0])]
        source, text = hit_tag(calls, model, hits["giants"][protein], hits["descriptions"]) if hits else ("", "")
        items.append([key[protein], model[0], model[1], round(score, 3), calls, source, text, []])
        taken.append(words)
    rule = ("2,000-aa models by Talea window score; sequential windows collapse to the higher-scoring one, "
            "near-identical copies in other giants fold into the first")
    payload = {"rule": rule + ("; listed spans and their copies left out" if excluded else ""), "items": items}
    (giants / "ranked.json").write_text(json.dumps(payload, separators=(",", ":")))
    print(f"{len(items)} entries ({sum(item[3] >= 1.5 for item in items)} at window score >= 1.5) from "
          f"{sum(map(len, scores.values()))} windows; {sum(len(item[7]) for item in items)} copies folded in; "
          f"left out: {left_out['excluded']} on listed spans, {left_out['copies of excluded']} copies of them, "
          f"{left_out['no model']} without a 2,000-aa model; entries with no repeat region: "
          f"{sum(not item[4] for item in items)}")


if __name__ == "__main__":
    main()
