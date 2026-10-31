"""Compact backbones for the giant-protein models: 5 bytes per residue, rebuilt in the browser.

A model is stored as its C-alpha trace on a 0.1 A grid (the first C-alpha as three int16, then int8 steps), one
byte per peptide giving the rotation of the peptide plane about the C-alpha(i) -> C-alpha(i+1) axis, and the C-alpha
pLDDT as a byte. Decoding places C(i), O(i) and N(i+1) of every peptide at the mean in-plane positions measured on
the models themselves (PLANE), which keeps the cartoon and the secondary structure 3Dmol assigns from N-O hydrogen
bonds. The page decodes with the same arithmetic (reference frame from the quantized trace alone).

Record layout (little-endian): uint32 n; int16[3] first C-alpha; uint16 k, then k x (uint32 step index, int16[3]
step) for steps longer than an int8 holds (12.7 A, model chain breaks); int8[3 (n-1)] steps (0 where listed above);
uint8[n-1] plane angles; uint8[n] pLDDT.
"""
import struct

import numpy as np

GRID = 10.0  # grid points per angstrom
# Mean in-plane coordinates (angstrom) of C(i), O(i), N(i+1) in the frame of C-alpha(i): along the C-alpha axis,
# then along the plane direction that points to O(i). Measured by measure_plane() on w500/w2000 models.
PLANE = {"C": (1.427, 0.523), "O": (1.656, 1.735), "N": (2.347, -0.372)}


def read_backbone(path):
    """N, CA, C, O coordinates, residue names and C-alpha B-factors (pLDDT) of a single-chain mmCIF model."""
    columns, atoms, index = [], {}, None
    with open(path) as handle:
        for line in handle:
            if line.startswith("_atom_site."):
                columns.append(line.strip().split(".", 1)[1])
            elif line.startswith("ATOM"):
                if index is None:
                    index = [columns.index(key) for key in (
                        "label_atom_id", "label_comp_id", "label_seq_id", "B_iso_or_equiv", "Cartn_x")]
                name_at, resn_at, seq_at, b_at, x_at = index
                fields = line.split()
                name = fields[name_at]
                if name in ("N", "CA", "C", "O"):
                    seq = int(fields[seq_at])
                    entry = atoms.setdefault(seq, {"resn": fields[resn_at]})
                    entry[name] = (float(fields[x_at]), float(fields[x_at + 1]), float(fields[x_at + 2]))
                    if name == "CA":
                        entry["b"] = float(fields[b_at])
    order = sorted(atoms)
    if order != list(range(order[0], order[0] + len(order))):
        raise ValueError(f"{path}: residue numbering has gaps")
    missing = [seq for seq in order if len(atoms[seq]) != 6]
    if missing:
        raise ValueError(f"{path}: incomplete backbone at {missing[:5]}")
    stack = lambda key: np.array([atoms[seq][key] for seq in order], dtype=np.float64)
    return {
        "N": stack("N"), "CA": stack("CA"), "C": stack("C"), "O": stack("O"),
        "resn": [atoms[seq]["resn"] for seq in order],
        "plddt": np.array([atoms[seq]["b"] for seq in order]),
    }



def frames(trace):
    """Unit C-alpha axis and two perpendicular reference vectors for every peptide, from the trace alone."""
    axis = np.diff(trace, axis=0)
    axis /= np.linalg.norm(axis, axis=1, keepdims=True)
    reference = np.zeros_like(axis)
    reference[np.arange(len(axis)), np.argmin(np.abs(axis), axis=1)] = 1.0
    first = reference - np.sum(reference * axis, axis=1, keepdims=True) * axis
    first /= np.linalg.norm(first, axis=1, keepdims=True)
    second = np.cross(axis, first)
    return axis, first, second


def encode(model):
    grid = np.rint(model["CA"] * GRID).astype(np.int64)
    steps = np.diff(grid, axis=0)
    if np.abs(grid).max() > 32767:
        raise ValueError("trace outside the int16 range")
    long_steps = np.flatnonzero(np.abs(steps).max(axis=1) > 127)
    if len(long_steps) > 65535:
        raise ValueError("too many long steps")
    trace = grid / GRID
    axis, first, second = frames(trace)
    toward_o = model["O"][:-1] - trace[:-1]
    toward_o -= np.sum(toward_o * axis, axis=1, keepdims=True) * axis
    angle = np.arctan2(np.sum(toward_o * second, axis=1), np.sum(toward_o * first, axis=1))
    plane = np.rint(angle / (2 * np.pi) * 256).astype(np.int64) % 256
    plddt = np.clip(np.rint(model["plddt"]), 0, 100).astype(np.uint8)
    n = len(grid)
    short = steps.copy()
    short[long_steps] = 0
    return b"".join([
        struct.pack("<I", n),
        grid[0].astype("<i2").tobytes(),
        struct.pack("<H", len(long_steps)),
        b"".join(struct.pack("<I", index) + steps[index].astype("<i2").tobytes() for index in long_steps),
        short.astype("<i1").tobytes(),
        plane.astype(np.uint8).tobytes(),
        plddt.tobytes(),
    ])


def decode(blob, offset=0):
    n = struct.unpack_from("<I", blob, offset)[0]
    position = offset + 4
    origin = np.frombuffer(blob, "<i2", 3, position).astype(np.int64)
    position += 6
    long_count = struct.unpack_from("<H", blob, position)[0]
    position += 2
    long_steps = []
    for _ in range(long_count):
        index = struct.unpack_from("<I", blob, position)[0]
        long_steps.append((index, np.frombuffer(blob, "<i2", 3, position + 4).astype(np.int64)))
        position += 10
    steps = np.frombuffer(blob, "<i1", 3 * (n - 1), position).astype(np.int64).reshape(-1, 3)
    for index, step in long_steps:
        steps[index] = step
    position += 3 * (n - 1)
    plane = np.frombuffer(blob, np.uint8, n - 1, position).astype(np.float64)
    position += n - 1
    plddt = np.frombuffer(blob, np.uint8, n, position).astype(np.float64)
    position += n
    trace = np.vstack([origin, origin + np.cumsum(steps, axis=0)]) / GRID
    axis, first, second = frames(trace)
    angle = plane / 256 * 2 * np.pi
    direction = np.cos(angle)[:, None] * first + np.sin(angle)[:, None] * second
    place = lambda atom: trace[:-1] + PLANE[atom][0] * axis + PLANE[atom][1] * direction
    carbon, oxygen, nitrogen = place("C"), place("O"), place("N")
    # Termini: N(0) and C/O(n-1) mirror their neighbours' geometry.
    n_first = trace[0] - (nitrogen[0] - trace[1]) if n > 1 else trace[0]
    c_last = trace[-1] + (trace[-2] - nitrogen[-1]) if n > 1 else trace[-1]
    o_last = c_last + (oxygen[-1] - carbon[-1])
    return {
        "CA": trace,
        "N": np.vstack([n_first, nitrogen]),
        "C": np.vstack([carbon, c_last]),
        "O": np.vstack([oxygen, o_last]),
        "plddt": plddt,
    }, position


def measure_plane(models):
    """Mean in-plane coordinates of C(i), O(i), N(i+1) over a set of models (for PLANE)."""
    sums = {atom: [] for atom in ("C", "O", "N")}
    for model in models:
        trace = model["CA"]
        axis, _, _ = frames(trace)
        toward_o = model["O"][:-1] - trace[:-1]
        toward_o -= np.sum(toward_o * axis, axis=1, keepdims=True) * axis
        direction = toward_o / np.linalg.norm(toward_o, axis=1, keepdims=True)
        for atom, coordinates in (("C", model["C"][:-1]), ("O", model["O"][:-1]), ("N", model["N"][1:])):
            offset = coordinates - trace[:-1]
            sums[atom].append(np.stack([np.sum(offset * axis, 1), np.sum(offset * direction, 1)], 1))
    return {atom: tuple(np.round(np.median(np.vstack(values), axis=0), 3)) for atom, values in sums.items()}
