# solreview

The Solenoid Review Bench is a browser page for checking predicted protein structures against
[Talea](https://github.com/jwestrob/Talea)'s repeat calls. For each model you record a solenoid call (yes, no or
unsure), a fold class, the repeat residues and notes.

- **Review sets** show one model per protein beside Talea's layer-13 and layer-35 attention maps. Blind mode keeps
  Talea's regions hidden until you have made the call. Where the dataset carries them, the page also shows Talea's
  score with its whole-protein and residue model z values, the predicted type, and annotation hits (Pfam, KOfam and
  others) for the model on screen. A set can be listed by Talea score or by the whole-protein minus residue z gap.
  Talea 0.3.2 adds regions from its residue model alone, where the attention search found no period. The page draws
  them as dotted outlines on the maps and in one flat colour on the structure. Each protein names the Talea version
  behind its regions.
- **Giants** shows long proteins that Talea scored in 2,000-aa windows. Each protein has:
  - a whole-protein track with Talea's windows and regions, the attention band and the residue state;
  - the folded models at every window size;
  - the Pfam, CDD and Foldseek hits for the model on screen;
  - a list of 2,000-aa models ranked by Talea window score, highest first.
- **Flagged** lists your own labelled proteins from every review set, newest label first, with each note under its
  row. A protein is flagged when its note contains the word *flag* ("FLAG!", "flag for refolding"). Two more buttons
  list every protein with a note, or everything you have labelled; labelled giant windows follow the list and open on
  the Giants tab. A call or note changed there saves to the protein's own set.

## Run it

```sh
git clone https://github.com/jwestrob/solreview.git
cd solreview
# put a review dataset in data/ (docs/data-format.md)
python3 serve.py --open
```

`serve.py` needs Python 3.9 or newer and its standard library. It serves the page at http://127.0.0.1:8765. The
page loads 3Dmol.js and its fonts from public CDNs.

Options:
- `--data DIR` reads a dataset from another folder.
- `--reviewer NAME` or `--labels FILE` chooses the labels file.
- `--port N` changes the port.

## Labels

Every change saves at once to `data/labels/<reviewer>.json`. The reviewer defaults to your user name, and the
status line under Notes names the file. Labels sit beside the dataset they describe.

The other files in `data/labels/` hold other reviewers' labels, which the page shows read-only:

- Where you have no call yet, a dashed pill in the list shows theirs.
- The sidebar shows each reviewer's call, fold and notes for the model on screen.
- With **Blind until I label** on, a review-set model shows them only after you make your own call.

To share labels, copy the files between `data/labels/` folders. Each reviewer writes only to their own file.
**Copy labels (TSV)** copies your labels to the clipboard as a table.

The file maps label ids to labels:

```json
{"census-review_01": {"solenoid": "yes", "fold": "alpha-solenoid", "residues": "120-340", "notes": "",
                      "talea_visible": false, "set": "census", "file": "review_01.pdb",
                      "updatedAt": "2026-10-30T21:53:47.387Z"},
 "giants-65-w2000_s1000-0": {"solenoid": "yes", "protein": "<record id>", "scheme": "w2000_s1000",
                             "start": 1, "end": 2000, "set": "giants", "updatedAt": "..."}}
```

- **Ids.** A review-set id is `<set>-<file stem>`. A giant id is `giants-<bundle key>-<scheme>-<model start>`, with
  the start 0-based.
- **Giant spans.** A giant label stores the model's span 1-based and inclusive.
- **`fold`** is one of `alpha-solenoid`, `beta-solenoid`, `alpha-beta-solenoid`, `helical-bundle`, `coiled-coil`,
  `other-fold` or `low-confidence`.
- **`talea_visible`** records whether Talea's regions were on screen when you made the solenoid call.

The same `index.html` also runs as a claude.ai artifact, where labels go to the artifact's database. Served by a plain
static file server, the page keeps labels in the browser.

## Data

Git ignores `data/`, labels included. [docs/data-format.md](docs/data-format.md) describes every file the page reads.
[build/README.md](build/README.md) explains how to build a dataset from folded models and Talea.
