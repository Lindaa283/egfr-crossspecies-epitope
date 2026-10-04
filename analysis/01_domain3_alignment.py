"""
Human vs mouse EGFR domain III comparison.

Domain III, mature numbering 310-481 (172 residues).

The mouse file is a precursor and carries a signal peptide; the human file
is already mature. Rather than hard-coding the offset, the mature start is
located by motif so the script fails loudly if a file is swapped.

Usage:  python analysis/01_domain3_alignment.py
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SEQ = ROOT / "data" / "sequences"

MATURE_START_MOTIF = "LEEKKVCQGTSN"
D3_FIRST, D3_LAST = 310, 481
D3_LEN = D3_LAST - D3_FIRST + 1


def load_mature(path):
    seq = "".join(path.read_text().split()).upper()
    offset = seq.find(MATURE_START_MOTIF)
    if offset < 0:
        raise SystemExit(f"{path.name}: mature start motif not found.")
    return seq[offset:]


def main():
    human = load_mature(SEQ / "human_EGFR_mature.txt")
    mouse = load_mature(SEQ / "mouse_EGFR_precursor.txt")

    human_d3 = human[D3_FIRST - 1:D3_LAST]
    mouse_d3 = mouse[D3_FIRST - 1:D3_LAST]

    for name, seg in (("human", human_d3), ("mouse", mouse_d3)):
        if len(seg) != D3_LEN:
            raise SystemExit(
                f"{name} domain III is {len(seg)} aa, expected {D3_LEN}. "
                "Check that the input file is the full sequence."
            )
        print(f"{name:>5} domain III: {seg[:10]} ... {seg[-10:]}  "
              f"({len(seg)} aa)")

    matches = sum(a == b for a, b in zip(human_d3, mouse_d3))
    print()
    print(f"identical positions : {matches}/{D3_LEN}")
    print(f"identity            : {matches / D3_LEN * 100:.1f}%")
    print(f"differences         : {D3_LEN - matches}")
    print()
    print("Differences (mature numbering, human>mouse):")

    diffs = []
    for position, (h, m) in enumerate(zip(human_d3, mouse_d3), start=D3_FIRST):
        if h != m:
            diffs.append((position, h, m))
            print(f"  {position}: {h}>{m}")

    # Where do the differences sit? Report density by 20-residue window so
    # the clustering claim is a measurement, not an impression.
    print()
    print("Difference density by window:")
    for start in range(D3_FIRST, D3_LAST + 1, 20):
        end = min(start + 19, D3_LAST)
        n = sum(1 for p, _, _ in diffs if start <= p <= end)
        bar = "#" * n
        print(f"  {start:>3}-{end:<3} {n:>2}  {bar}")


if __name__ == "__main__":
    main()
