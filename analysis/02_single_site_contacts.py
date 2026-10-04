"""
Contacts between one EGFR chain A residue and the antibody chains.

PDB file is not committed (1.3 MB). Download it first:
    curl -o data/structures/6ARU.pdb https://files.rcsb.org/download/6ARU.pdb

Note: 6ARU chains B and C are an antibody. Contacts measured here are used
to calibrate which EGFR surface positions are interface-capable. No part of
that antibody is used in any design.

Usage:  python analysis/02_single_site_contacts.py
"""

from pathlib import Path
from math import dist

ROOT = Path(__file__).resolve().parent.parent
pdb = ROOT / "data" / "structures" / "6ARU.pdb"
if not pdb.exists():
    raise SystemExit(f"Missing {pdb} - see the download command above.")

atoms = []

for line in pdb.read_text(errors="ignore").splitlines():
    if not line.startswith("ATOM"):
        continue

    chain = line[21].strip()
    resnum = line[22:26].strip()
    resname = line[17:20].strip()
    element = line[76:78].strip().upper()

    # Exclude hydrogens; compare protein heavy atoms only
    if element in {"H", "D"}:
        continue

    xyz = (
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54]),
    )
    atoms.append((chain, resnum, resname, xyz))

# Which chain A residue to inspect (mature numbering).
SITE = "471"

selected = [a for a in atoms if a[0] == "A" and a[1] == SITE]
if not selected:
    raise SystemExit("Residue not found in chain A - check PDB numbering.")

for target_chain in ("B", "C"):
    target = [a for a in atoms if a[0] == target_chain]
    if not target:
        print(f"没找到 Chain {target_chain}")
        continue

    contacts = {}
    for a in selected:
        for b in target:
            d = dist(a[3], b[3])
            key = (b[1], b[2])
            contacts[key] = min(contacts.get(key, float("inf")), d)

    nearest = min(contacts.items(), key=lambda item: item[1])
    close = [(key, d) for key, d in contacts.items() if d < 5.0]

    print(
        f"到链 {target_chain} 最近：残基 {nearest[0][1]} "
        f"{nearest[0][0]}，{nearest[1]:.2f} Å；"
        f"5 Å 内有 {len(close)} 个抗体残基"
    )