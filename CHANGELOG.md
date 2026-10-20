# Changelog

## 1.0.0 (2026-10-20)

First release of the combined model used in Pettebone, Han et al., bioRxiv 2026
(DOI 10.64898/2026.10.19.760310), Fig. 8A and Supplementary Fig. S5.

- hiMCES2026 (unchanged host) with ChR2 (Williams et al. 2013) and the modified
  Piezo1 model (two gates; activation and availability; three fitted kinetic
  parameters), as used for the preprint.
- Entry scripts: `MasterCompute_OptoSono.m` (Fig. 8A and controls),
  `FigureS5_Piezo1.m` (Fig. S5), `condition_OptoSono.m` (conditioning),
  `verify_release.m` (checks and reproduction).
- Conditioned states, conditioning start states, reference traces and reference
  channel measurements of the preprint.
- Optional calibration scripts (`calibration/`), which write only to `results/`.

Relation to the development code: equations, parameters, initial conditions,
solver settings and protocols are those used for the preprint. Files were renamed
for publication (development name of the modified Piezo1 model: "variant B"):
`hiMCES2026_piezo1B_chr2.m` to `model/hiMCES2026_ChR2_Piezo1.m`,
`piezo1B_chr2_opts.m` to `model/ChR2_Piezo1_opts.m`, `piezo1B_channel.m` to
`model/Piezo1_channel.m`, and the calibration functions `targets_B`, `measure_B`,
`fit_B`, `check_B` to `Piezo1_targets`, `measure_Piezo1`, `fit_Piezo1`,
`check_Piezo1`. In the model file only comments and the function name changed.
The optional fit now includes the restart step that produced the frozen values and
writes to `results/`. `verify_release.m` reproduces the stored preprint runs with
zero difference in the final states on the computer that produced them.
