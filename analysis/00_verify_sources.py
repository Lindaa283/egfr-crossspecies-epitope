"""
Source verification. Run this first.

Checks three things before any downstream analysis is trusted:
  1. The mouse file is a precursor (includes signal peptide); find the offset
     to mature numbering rather than assuming it.
  2. The human file and the 6ARU chain A construct are not identical;
     report exactly where they differ and whether those positions fall
     inside the region we analyse.
  3. The domain III slices really are 172 residues and really start at
     mature position 310 in both species.

Usage:  python analysis/00_verify_sources.py
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SEQ = ROOT / "data" / "sequences"
STRUCT = ROOT / "data" / "structures"

# Mature human EGFR begins with this motif. Used to locate the mature
# start in a precursor sequence instead of hard-coding an offset.
MATURE_START_MOTIF = "LEEKKVCQGTSN"

# Domain III, mature numbering (inclusive).
D3_FIRST, D3_LAST = 310, 481


def read_plain(path):
    return "".join(path.read_text().split()).upper()


def read_fasta(path):
    return "".join(
        line.strip()
        for line in path.read_text().splitlines()
        if not line.startswith(">")
    ).upper()


def main():
    human = read_plain(SEQ / "human_EGFR_mature.txt")
    mouse = read_plain(SEQ / "mouse_EGFR_precursor.txt")
    pdb_seq = read_fasta(STRUCT / "6aru_chain_A.fasta")

    print("=" * 62)
    print("1. Lengths")
    print("=" * 62)
    print(f"  human file        : {len(human)} aa")
    print(f"  mouse file        : {len(mouse)} aa")
    print(f"  6ARU chain A      : {len(pdb_seq)} aa")

    print()
    print("=" * 62)
    print("2. Mature-sequence offsets")
    print("=" * 62)
    human_offset = human.find(MATURE_START_MOTIF)
    mouse_offset = mouse.find(MATURE_START_MOTIF)
    for name, offset in (("human", human_offset), ("mouse", mouse_offset)):
        if offset < 0:
            print(f"  {name}: mature start motif NOT FOUND — check the file")
            return
        print(f"  {name}: mature residue 1 is at index {offset}", end="")
        print(" (signal peptide present)" if offset else " (already mature)")

    print()
    print("=" * 62)
    print("3. human file vs 6ARU chain A construct")
    print("=" * 62)
    overlap = min(len(human), len(pdb_seq))
    diffs = [
        (i + 1, human[i], pdb_seq[i])
        for i in range(overlap)
        if human[i] != pdb_seq[i]
    ]
    if not diffs:
        print("  identical across the overlap")
    else:
        print(f"  {len(diffs)} mismatch(es) in the first {overlap} positions:")
        for pos, h, p in diffs:
            inside = D3_FIRST <= pos <= D3_LAST
            flag = "  <-- INSIDE domain III" if inside else ""
            print(f"    mature {pos:>4}:  file={h}  6ARU={p}{flag}")
    if len(pdb_seq) > overlap:
        print(f"  6ARU has {len(pdb_seq) - overlap} extra residue(s): "
              f"{pdb_seq[overlap:]}  (purification tag, not native)")

    affected = [p for p, _, _ in diffs if D3_FIRST <= p <= D3_LAST]
    print()
    if affected:
        print(f"  WARNING: {len(affected)} difference(s) fall inside "
              f"domain III — the alignment result is affected.")
    else:
        print("  All differences fall OUTSIDE domain III "
              f"({D3_FIRST}-{D3_LAST}); the alignment is unaffected.")

    print()
    print("=" * 62)
    print("4. Domain III slice check")
    print("=" * 62)
    expected_len = D3_LAST - D3_FIRST + 1
    h_start = human_offset + D3_FIRST - 1
    m_start = mouse_offset + D3_FIRST - 1
    h_d3 = human[h_start:h_start + expected_len]
    m_d3 = mouse[m_start:m_start + expected_len]

    ok = True
    for name, seg in (("human", h_d3), ("mouse", m_d3)):
        status = "OK" if len(seg) == expected_len else "LENGTH MISMATCH"
        if len(seg) != expected_len:
            ok = False
        print(f"  {name}: {len(seg)} aa [{status}]  "
              f"{seg[:10]}...{seg[-10:]}")

    if not ok:
        print("  Stop: a slice is the wrong length. Check the input files.")
        return

    identical = sum(a == b for a, b in zip(h_d3, m_d3))
    print()
    print(f"  identity: {identical}/{expected_len} "
          f"= {identical / expected_len * 100:.1f}%")
    print(f"  differences: {expected_len - identical}")

    print()
    print("All source checks completed.")


if __name__ == "__main__":
    main()
