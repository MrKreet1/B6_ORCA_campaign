# B6 vibrational analysis

This folder contains vibrational analysis results for the selected best B6 structure.

## Selected structure

- Cluster: B6
- Charge: 0
- Multiplicity: 3
- Method: PBE0-D4/def2-TZVP OptFreq
- Imaginary frequencies: 0
- Number of real vibrational modes: 12
- Frequency range: 233.78–1407.26 cm^-1

## Files

| File | Description |
|---|---|
| B6_all_vibrational_frequencies.csv | All 12 non-zero vibrational frequencies |
| B6_normal_mode_amplitudes.csv | Normal-mode displacement components dx, dy, dz and amplitude for each atom |
| B6_mode_summary.csv | Summary table: max amplitude and dominant atom for each mode |
| B6_vibrational_frequencies_raw.txt | Raw ORCA frequency block |
| B6_best.out | ORCA output file of the selected best structure |
| B6_best.hess | ORCA Hessian file containing vibrational information |
| B6_best_optimized.xyz | Optimized geometry of the selected B6 structure |

## Note

The displacement amplitudes are normalized normal-mode components from ORCA.
They describe relative participation of atoms in each vibrational mode and should not be interpreted as absolute thermal amplitudes in angstroms.