"""Write annotation hits for review-set models to BENCH/data/review_hits.json.

usage: python3 review_hits.py BENCH HITS_TSV [--sources "one line on where the hits come from"]

HITS_TSV has the columns item_id, source, start, end, family, description, evalue. item_id is the page's model id
(`<set>-<file stem>`). start and end are 0-based end-exclusive positions along the modelled sequence; leave both
empty for a whole-protein call (a KOfam ortholog, say). These are observed annotations, shown as they are: the page
names no family or system on its own.

Items in HITS_TSV replace their earlier entries; other items are kept.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("bench", type=Path)
    parser.add_argument("hits", type=Path)
    parser.add_argument("--sources", default="")
    args = parser.parse_args()
    out = args.bench / "data" / "review_hits.json"
    payload = json.loads(out.read_text()) if out.exists() else {"sources": "", "items": {}}
    fresh: dict[str, list] = {}
    with args.hits.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            start = int(row["start"]) if row["start"] else None
            end = int(row["end"]) if row["end"] else None
            evalue = float(row["evalue"]) if row["evalue"] else None
            fresh.setdefault(row["item_id"], []).append(
                [row["source"], start, end, row["family"], row["description"], evalue])
    payload["items"].update(fresh)
    if args.sources:
        payload["sources"] = args.sources
    out.write_text(json.dumps(payload, separators=(",", ":")))
    print(f"{len(fresh)} items, {sum(len(v) for v in fresh.values())} hits; {len(payload['items'])} items in {out}")


if __name__ == "__main__":
    main()
