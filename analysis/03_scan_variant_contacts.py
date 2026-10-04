"""
Scan all 21 human/mouse domain III variant positions and measure how
close each sits to the antibody interface in 6ARU.

PDB file is not committed (1.3 MB). Download it first:
    curl -o data/structures/6ARU.pdb https://files.rcsb.org/download/6ARU.pdb

Note: 6ARU chains B and C are an antibody. Contacts measured here are used
to calibrate which EGFR surface positions are interface-capable. No part of
that antibody is used in any design.

Usage:  python analysis/03_scan_variant_contacts.py
"""

from pathlib import Path
from math import dist

ROOT = Path(__file__).resolve().parent.parent
pdb = ROOT / "data" / "structures" / "6ARU.pdb"
if not pdb.exists():
    raise SystemExit(f"Missing {pdb} - see the download command above.")


# 成熟 EGFR 编号：人源氨基酸 > 小鼠氨基酸
variants = {
    324: "S>T", 337: "N>Y", 340: "S>A", 353: "R>K",
    359: "H>R", 366: "Q>R", 369: "D>E", 388: "E>D",
    390: "R>W", 418: "S>G", 443: "K>R", 460: "S>P",
    461: "G>N", 467: "I>M", 468: "S>N", 471: "G>A",
    473: "N>K", 474: "S>D", 478: "T>V", 479: "G>N",
    480: "Q>H",
}

aa3_to_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

atoms = []
for line in pdb.read_text(errors="ignore").splitlines():
    if not line.startswith("ATOM"):
        continue

    element = line[76:78].strip().upper()
    if element in {"H", "D"}:
        continue

    atoms.append({
        "chain": line[21].strip(),
        "resnum": line[22:26].strip(),
        "resname": line[17:20].strip(),
        "xyz": (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        ),
    })

def measure(site_atoms, chain_atoms):
    residue_min = {}

    for a in site_atoms:
        for b in chain_atoms:
            key = (b["resnum"], b["resname"])
            d = dist(a["xyz"], b["xyz"])
            residue_min[key] = min(residue_min.get(key, float("inf")), d)

    if not residue_min:
        return "no atoms", 0

    nearest = min(residue_min.items(), key=lambda item: item[1])
    nearest_text = f'{nearest[0][1]}{nearest[0][0]} {nearest[1]:.2f} A'
    within_5 = sum(d < 5.0 for d in residue_min.values())
    return nearest_text, within_5

for position, change in variants.items():
    site_atoms = [
        a for a in atoms
        if a["chain"] == "A" and a["resnum"] == str(position)
    ]

    if not site_atoms:
        print(f"{position} {change}: 没找到结构残基")
        continue

    observed = aa3_to_1.get(site_atoms[0]["resname"], "?")
    expected = change.split(">")[0]
    check = "OK" if observed == expected else f"核对！结构里是 {observed}"

    b_result = measure(
        site_atoms, [a for a in atoms if a["chain"] == "B"]
    )
    c_result = measure(
        site_atoms, [a for a in atoms if a["chain"] == "C"]
    )

    print(
        f"{position} {change} [{check}] | "
        f"B: 最近 {b_result[0]}, <5 A 残基数 {b_result[1]} | "
        f"C: 最近 {c_result[0]}, <5 A 残基数 {c_result[1]}"
    )