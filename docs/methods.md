# Methods

Written as work is actually performed. Nothing is recorded here before it runs.

## 1. Sources and retrieval

| Item | Accession | Form |
|---|---|---|
| Human EGFR | UniProt P00533 | mature ectodomain, 621 aa |
| Mouse EGFR | UniProt Q01279 | full precursor, 1210 aa |
| Structure | PDB 6ARU | chain A = EGFR construct; chains B, C = antibody |

Retrieval dates and checksums: *(to record)*

## 2. Numbering convention

All positions are **mature human EGFR numbering**.

The mouse record is a precursor with a 24-residue signal peptide. Mature
residue 1 is located by the motif `LEEKKVCQGTSN` rather than by a hard-coded
offset, so a substituted or re-downloaded file cannot silently shift the
frame.

Domain III is taken as mature 310-481 (172 residues).

## 3. Construct vs native sequence

6ARU chain A is a crystallisation construct. Against the reference human
sequence it differs at mature 516 (N->K) and 610 (E->R), and carries a
C-terminal polyhistidine tag. All differences lie outside domain III, so
the cross-species comparison is unaffected. The check is run by
`analysis/00_verify_sources.py` on every execution rather than assumed.

## 4. Cross-species comparison

Direct position-by-position comparison of the two 172-residue domain III
segments. No gap-opening alignment was needed: the segments are equal length
and ungapped across this region.

Result: 151/172 identical (87.8%), 21 differences.

Difference positions are reported by the script rather than transcribed by
hand. They were independently cross-checked against the hand-entered list
used in `analysis/03_scan_variant_contacts.py`; the two agree exactly.

## 5. Structural context of the variant positions

Heavy-atom minimum distances from each variant position in chain A to the
antibody chains B and C in 6ARU. Hydrogens and deuteriums excluded. A
residue is counted as within contact range at < 5.0 A.

Purpose: calibration only — to establish which EGFR surface positions are
capable of sitting at a protein-protein interface. The antibody itself is
not used as a design starting point.

Each position is validated before measurement: the three-letter residue name
present in the structure is converted to one-letter and compared against the
residue the sequence alignment predicts. A mismatch is reported rather than
measured.

## 6. Epitope definition and cross-species conservation

The structural epitope is defined as every chain A residue with at least one
heavy atom within 5.0 A of chains B or C. This yields 26 residues, all within
domain III.

Intersecting that set with the 21 species differences: 19 conserved, 7
divergent (27%). Four of the seven fall within residues 465-473. The two
shortest contacts in the entire epitope are both species differences
(S468N, 2.82 A; I467M, 3.05 A).

Validation: cetuximab is experimentally documented as human-specific and not
cross-reactive with mouse EGFR. This analysis uses no binding data and
recovers a structural basis consistent with that observation, which is
treated here as a sanity check on the method rather than as a new finding.

Scope: 6ARU describes one antibody's binding surface. The result establishes
that these positions can participate in a protein-protein interface; it does
not establish that any arbitrary binder must contact them.

## 6b. Epitope selection for design

*(to write - the reasoning from the data above to a chosen target surface)*

## 7. Design generation

*(to write — model, version, parameters, number of candidates generated,
number filtered out, and why)*

## 8. Limitations

- Conservation at the sequence level does not guarantee conservation of
  surface conformation; identical residues can present differently.
- Contact distances from a single crystal structure describe one bound
  state, not the accessible conformational ensemble.
- Nothing in this repository is experimental evidence. No computational
  score establishes binding, affinity, cross-reactivity, or pH selectivity.
