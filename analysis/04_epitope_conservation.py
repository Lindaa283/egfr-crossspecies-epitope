"""
Define the cetuximab epitope on EGFR from 6ARU, then ask how much of it is
conserved between human and mouse.

6ARU is a cetuximab Fab mutant bound to the EGFR extracellular region.
Chain A is EGFR; chains B and C are the Fab. This script:

  1. Collects every chain A residue with a heavy atom within CUTOFF of the
     Fab -- that set is the structural epitope.
  2. Intersects it with the 21 human/mouse domain III differences.
  3. Reports the conserved fraction of the epitope.

Calibration only. Cetuximab defines ONE binding surface on EGFR, not the
only one. No part of this antibody is used in any design.

Usage:  python analysis/04_epitope_conservation.py
"""

from pathlib import Path
from math import dist

ROOT = Path(__file__).resolve().parent.parent
PDB = ROOT / "data" / "structures" / "6ARU.pdb"
SEQ = ROOT / "data" / "sequences"

CUTOFF = 5.0
EGFR_CHAIN = "A"
FAB_CHAINS = ("B", "C")

MATURE_START_MOTIF = "LEEKKVCQGTSN"
D3_FIRST, D3_LAST = 310, 481

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def load_mature(path):
    seq = "".join(path.read_text().split()).upper()
    offset = seq.find(MATURE_START_MOTIF)
    if offset < 0:
        raise SystemExit(f"{path.name}: mature start motif not found.")
    return seq[offset:]


def load_atoms():
    if not PDB.exists():
        raise SystemExit(
            f"Missing {PDB}\n"
            "  curl -o data/structures/6ARU.pdb "
            "https://files.rcsb.org/download/6ARU.pdb"
        )
    atoms = []
    for line in PDB.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[76:78].strip().upper() in {"H", "D"}:
            continue
        atoms.append({
            "chain": line[21].strip(),
            "resnum": line[22:26].strip(),
            "resname": line[17:20].strip(),
            "xyz": (float(line[30:38]), float(line[38:46]),
                    float(line[46:54])),
        })
    return atoms


def main():
    human = load_mature(SEQ / "human_EGFR_mature.txt")
    mouse = load_mature(SEQ / "mouse_EGFR_precursor.txt")

    human_d3 = human[D3_FIRST - 1:D3_LAST]
    mouse_d3 = mouse[D3_FIRST - 1:D3_LAST]
    variants = {
        pos: (h, m)
        for pos, (h, m) in enumerate(zip(human_d3, mouse_d3), start=D3_FIRST)
        if h != m
    }

    atoms = load_atoms()
    egfr = [a for a in atoms if a["chain"] == EGFR_CHAIN]
    fab = [a for a in atoms if a["chain"] in FAB_CHAINS]
    if not egfr or not fab:
        raise SystemExit("Expected chains not found in the structure.")

    # Minimum heavy-atom distance from each EGFR residue to the Fab.
    nearest = {}
    for a in egfr:
        key = (int(a["resnum"]), a["resname"])
        best = nearest.get(key, float("inf"))
        for b in fab:
            d = dist(a["xyz"], b["xyz"])
            if d < best:
                best = d
        nearest[key] = best

    epitope = sorted(
        (pos, name, d) for (pos, name), d in nearest.items() if d < CUTOFF
    )

    print("=" * 66)
    print(f"Cetuximab epitope on EGFR  (chain A residues within "
          f"{CUTOFF} A of the Fab)")
    print("=" * 66)
    print(f"  epitope size: {len(epitope)} residues")
    in_d3 = [e for e in epitope if D3_FIRST <= e[0] <= D3_LAST]
    print(f"  of which in domain III ({D3_FIRST}-{D3_LAST}): {len(in_d3)}")
    print()

    print(f"{'pos':>5}  {'aa':<4} {'min dist':>9}   species difference")
    print("-" * 66)
    conserved = divergent = 0
    for pos, name, d in epitope:
        one = AA3.get(name, "?")
        if pos in variants:
            h, m = variants[pos]
            note = f"DIFFERS  {h}>{m}"
            divergent += 1
            if one != h:
                note += f"   [CHECK: structure has {one}]"
        else:
            note = "conserved"
            conserved += 1
        print(f"{pos:>5}  {one}{name:<3} {d:>7.2f} A   {note}")

    print("-" * 66)
    total = conserved + divergent
    print(f"  conserved: {conserved}/{total} "
          f"({conserved / total * 100:.0f}%)")
    print(f"  divergent: {divergent}/{total} "
          f"({divergent / total * 100:.0f}%)")

    print()
    print("=" * 66)
    print("The 21 domain III differences, ranked by closest approach")
    print("=" * 66)
    rows = []
    for pos, (h, m) in variants.items():
        key = next((k for k in nearest if k[0] == pos), None)
        d = nearest.get(key, float("inf")) if key else float("inf")
        rows.append((d, pos, h, m))
    rows.sort()
    for d, pos, h, m in rows:
        mark = "AT INTERFACE" if d < CUTOFF else ""
        shown = f"{d:>7.2f} A" if d != float("inf") else "   absent"
        print(f"  {pos:>4} {h}>{m}  {shown}   {mark}")

    at_interface = sum(1 for d, *_ in rows if d < CUTOFF)
    print()
    print(f"  {at_interface} of {len(rows)} species differences sit inside "
          f"this epitope.")


if __name__ == "__main__":
    main()
