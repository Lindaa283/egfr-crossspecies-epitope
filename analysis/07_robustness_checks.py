"""
Robustness checks for the two arbitrary parameters in this analysis.

Two numbers were chosen by convention rather than derived:

  1. The 5.0 A contact cutoff used to define the epitope (analysis/04).
  2. The domain III boundaries, 310-481 (analysis/01).

Neither is a physical constant. A result that only holds at one setting is an
artifact of that setting, so both are varied here and the conclusions are
re-derived at each value.

Usage:  python analysis/07_robustness_checks.py
"""

from pathlib import Path
from math import dist

ROOT = Path(__file__).resolve().parent.parent
PDB = ROOT / "data" / "structures" / "6ARU.pdb"
SEQ = ROOT / "data" / "sequences"

MATURE_START_MOTIF = "LEEKKVCQGTSN"
D3_FIRST, D3_LAST = 310, 481

CUTOFFS = [3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
BOUNDARIES = [
    (300, 490), (305, 485), (309, 481), (310, 481),
    (310, 480), (311, 480), (315, 475), (320, 470),
]


def load_mature(path):
    seq = "".join(path.read_text().split()).upper()
    i = seq.find(MATURE_START_MOTIF)
    if i < 0:
        raise SystemExit(f"{path.name}: mature start motif not found.")
    return seq[i:]


def variants(human, mouse, first, last):
    return {
        p: (h, m)
        for p, (h, m) in enumerate(
            zip(human[first - 1:last], mouse[first - 1:last]), start=first
        ) if h != m
    }


def nearest_to_fab():
    """Minimum heavy-atom distance from each chain A residue to the Fab."""
    if not PDB.exists():
        raise SystemExit(
            f"Missing {PDB}\n  curl -o data/structures/6ARU.pdb "
            "https://files.rcsb.org/download/6ARU.pdb"
        )
    egfr, fab = [], []
    for line in PDB.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[76:78].strip().upper() in {"H", "D"}:
            continue
        chain = line[21].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if chain == "A":
            egfr.append((int(line[22:26]), xyz))
        elif chain in ("B", "C"):
            fab.append(xyz)

    out = {}
    for num, xyz in egfr:
        best = out.get(num, float("inf"))
        for b in fab:
            d = dist(xyz, b)
            if d < best:
                best = d
        out[num] = best
    return out


def main():
    human = load_mature(SEQ / "human_EGFR_mature.txt")
    mouse = load_mature(SEQ / "mouse_EGFR_precursor.txt")

    # --- Check 1: contact cutoff -------------------------------------
    var = variants(human, mouse, D3_FIRST, D3_LAST)
    near = nearest_to_fab()

    print("=" * 72)
    print("Check 1 - does the epitope result depend on the 5 A cutoff?")
    print("=" * 72)
    print(f"{'cutoff':>7}  {'epitope':>8}  {'divergent':>10}  {'%':>5}   "
          "positions")
    print("-" * 72)
    always = None
    for cut in CUTOFFS:
        ep = {p for p, d in near.items() if d < cut}
        both = sorted(ep & set(var))
        always = set(both) if always is None else always & set(both)
        labels = " ".join(f"{var[p][0]}{p}{var[p][1]}" for p in both)
        print(f"{cut:>6.1f}A  {len(ep):>8}  {len(both):>10}  "
              f"{len(both) / len(ep) * 100:>4.0f}%   {labels}")
    print("-" * 72)
    print(f"  present at every cutoff: "
          f"{' '.join(f'{var[p][0]}{p}{var[p][1]}' for p in sorted(always))}")

    # --- Check 2: domain boundaries ----------------------------------
    print()
    print("=" * 72)
    print("Check 2 - does the identity figure depend on where domain III "
          "is cut?")
    print("=" * 72)
    print(f"{'range':>12}  {'length':>7}  {'identical':>10}  {'identity':>9}")
    print("-" * 72)
    values = []
    for first, last in BOUNDARIES:
        h = human[first - 1:last]
        m = mouse[first - 1:last]
        same = sum(a == b for a, b in zip(h, m))
        pct = same / len(h) * 100
        values.append(pct)
        print(f"{first}-{last:<7}  {len(h):>7}  {same:>10}  {pct:>8.1f}%")
    print("-" * 72)
    print(f"  range across all boundary definitions: "
          f"{min(values):.1f}% - {max(values):.1f}%")
    print()
    print("  Identity rises as the window is narrowed, which is consistent")
    print("  with the differences being concentrated at the two ends - the")
    print("  same pattern the window-density count in analysis/01 shows.")


if __name__ == "__main__":
    main()
