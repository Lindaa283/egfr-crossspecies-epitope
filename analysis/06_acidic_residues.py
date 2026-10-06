"""
Find solvent-exposed acidic residues (Asp, Glu) in the conserved stretch of
EGFR domain III.

Why: a pH-switchable binder needs histidines at the interface that are
protonated at pH 6.5 and neutral at pH 7.4. A protonated His carries positive
charge, so it needs a negatively charged partner on the target to pair with.
Asp and Glu carry that charge. If the conserved region has no exposed acidic
residues, the pH switch cannot be built there.

Method: solvent accessibility by Shrake-Rupley, implemented here rather than
imported, so the assumptions are visible. Each atom is given a van der Waals
radius, expanded by a 1.4 A probe (roughly a water molecule). Points are
spread evenly over that expanded sphere; a point counts as accessible if no
neighbouring atom's expanded sphere contains it. Residue SASA is the sum over
its atoms, and relative accessibility divides that by a reference maximum for
that residue type.

Chain A is measured alone, with the antibody removed, so what is reported is
the free receptor surface rather than the bound complex.

Usage:  python analysis/06_acidic_residues.py
"""

from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
PDB = ROOT / "data" / "structures" / "6ARU.pdb"
SEQ = ROOT / "data" / "sequences"

# The conserved stretch identified in analysis/01 (mature numbering).
REGION_FIRST, REGION_LAST = 370, 450
D3_FIRST, D3_LAST = 310, 481
MATURE_START_MOTIF = "LEEKKVCQGTSN"

PROBE = 1.4          # radius of a water molecule, Angstrom
N_POINTS = 256       # sphere sampling density
EXPOSED_CUTOFF = 0.25   # relative accessibility above this counts as surface

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}

# Reference maximum SASA per residue type (Tien et al. 2013, theoretical).
MAX_ACC = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
    "GLU": 223, "GLN": 225, "GLY": 104, "HIS": 224, "ILE": 197,
    "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
    "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174,
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def sphere_points(n):
    """Evenly spread n points on a unit sphere (golden spiral)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([
        np.cos(theta) * np.sin(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(phi),
    ])


def load_chain(path, chain):
    """Heavy atoms of one chain: (resnum, resname, element, xyz)."""
    atoms = []
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21].strip() != chain:
            continue
        element = line[76:78].strip().upper()
        if element in {"H", "D"} or element not in VDW:
            continue
        atoms.append((
            int(line[22:26]),
            line[17:20].strip(),
            element,
            (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        ))
    return atoms


def shrake_rupley(atoms):
    """Per-atom SASA in square Angstrom."""
    xyz = np.array([a[3] for a in atoms])
    radii = np.array([VDW[a[2]] for a in atoms]) + PROBE
    unit = sphere_points(N_POINTS)
    areas = np.zeros(len(atoms))

    # Only atoms within reach of the largest sphere pair can occlude.
    max_reach = radii.max() * 2
    for i in range(len(atoms)):
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        near = np.where((d < max_reach) & (d > 0))[0]
        test = xyz[i] + unit * radii[i]
        blocked = np.zeros(N_POINTS, dtype=bool)
        for j in near:
            blocked |= np.linalg.norm(test - xyz[j], axis=1) < radii[j]
        accessible = (~blocked).sum() / N_POINTS
        areas[i] = accessible * 4 * np.pi * radii[i] ** 2
    return areas


def load_mature(path):
    seq = "".join(path.read_text().split()).upper()
    i = seq.find(MATURE_START_MOTIF)
    if i < 0:
        raise SystemExit(f"{path.name}: mature start motif not found.")
    return seq[i:]


def main():
    if not PDB.exists():
        raise SystemExit(
            f"Missing {PDB}\n  curl -o data/structures/6ARU.pdb "
            "https://files.rcsb.org/download/6ARU.pdb"
        )

    human = load_mature(SEQ / "human_EGFR_mature.txt")
    mouse = load_mature(SEQ / "mouse_EGFR_precursor.txt")
    variants = {
        p for p, (h, m) in enumerate(
            zip(human[D3_FIRST - 1:D3_LAST], mouse[D3_FIRST - 1:D3_LAST]),
            start=D3_FIRST,
        ) if h != m
    }

    atoms = load_chain(PDB, "A")
    print(f"chain A heavy atoms: {len(atoms)}")
    print("computing solvent accessibility...")
    areas = shrake_rupley(atoms)

    # Sum per residue.
    res_area, res_name = {}, {}
    for (num, name, _, _), a in zip(atoms, areas):
        res_area[num] = res_area.get(num, 0.0) + a
        res_name[num] = name

    rows = []
    for num in sorted(res_area):
        if not (REGION_FIRST <= num <= REGION_LAST):
            continue
        name = res_name[num]
        if name not in ("ASP", "GLU"):
            continue
        rel = res_area[num] / MAX_ACC[name]
        rows.append((num, AA3[name], res_area[num], rel, num in variants))

    print()
    print("=" * 68)
    print(f"Acidic residues in the conserved stretch "
          f"{REGION_FIRST}-{REGION_LAST}")
    print("=" * 68)
    print(f"{'pos':>5}  {'aa':<3} {'SASA':>8}  {'rel':>6}   status")
    print("-" * 68)
    exposed = []
    for num, aa, area, rel, is_var in rows:
        if rel >= EXPOSED_CUTOFF:
            status = "EXPOSED"
            exposed.append((num, aa, rel))
        elif rel >= 0.10:
            status = "partly buried"
        else:
            status = "buried"
        if is_var:
            status += "   [differs in mouse]"
        print(f"{num:>5}  {aa:<3} {area:>7.1f} A2 {rel:>6.2f}   {status}")

    print("-" * 68)
    print(f"  acidic residues in region : {len(rows)}")
    print(f"  of which exposed (>={EXPOSED_CUTOFF:.0%}) : {len(exposed)}")

    if not exposed:
        print()
        print("  No exposed acidic residue here. A His-based pH switch would")
        print("  have no charged partner in this region.")
        return

    # How close are the exposed ones to each other? A His pair needs partners
    # within reach of one binding surface.
    print()
    print("=" * 68)
    print("Spatial clustering of the exposed acidic residues")
    print("=" * 68)
    coords = {}
    for line in PDB.read_text(errors="ignore").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" \
                and line[21].strip() == "A":
            coords[int(line[22:26])] = np.array([
                float(line[30:38]), float(line[38:46]), float(line[46:54])
            ])

    nums = [n for n, _, _ in exposed]
    print(f"{'':>7}" + "".join(f"{n:>7}" for n in nums))
    for a in nums:
        row = f"{a:>7}"
        for b in nums:
            if a == b:
                row += f"{'-':>7}"
            else:
                row += f"{np.linalg.norm(coords[a] - coords[b]):>7.1f}"
        print(row)
    print()
    print("  CA-CA distances in Angstrom. Two acidic residues within roughly")
    print("  15 A could be reached by one designed binding surface.")


if __name__ == "__main__":
    main()
