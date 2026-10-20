# Notice

This release (version 1.0.0, DOI 10.5281/zenodo.23430167) is distributed under the
MIT License in `LICENSE`, Copyright (c) 2026 Amin Forouzandehmehr. It contains or
builds on the components below. Their notices are retained.

## hiMCES2026 (host model)

- `upstream/hiMCES2026.m` is an unchanged copy of the published model file from
  https://github.com/aminforouzandehmehr/Forouzandehmehr2026-hiPSC-CMs-Model
  (commit 9a490f5b61a3d1f63bcf094b3e61948632f905ea; SHA-256 `70fb0c6b57c801378caebddd99c68ca880ee41338ed0c4e255167fc983a36306`).
  `model/heav.m` and `model/passiveForces.m` are unchanged copies from the same
  release. MIT License, Copyright (c) 2026 Amin Forouzandehmehr
  (`upstream/LICENSE_hiMCES2026`).
- Paper: Li W, Forouzandehmehr MA, McLeod D, Pozo MR, Heinson YW, Kay MW, Li Z,
  Morotti S, Entcheva E. Response of human iPSC-cardiomyocytes to adrenergic drugs
  assessed by high-throughput pericellular oxygen measurements and computational
  modeling. J Mol Cell Cardiol Plus 2026;15:100834. DOI 10.1016/j.jmccpl.2026.100834.
- The host equations and parameters inside `model/hiMCES2026_ChR2_Piezo1.m` are those
  of the published file. `verify_release.m` (check 1a) confirms that the combined
  model with the added currents off equals the published file in every derivative
  and output.

## Piezo1 (CinC 2023, code version 1.1)

- `upstream/Piezo1.m` is an unchanged copy of `Piezo1.m`, version 1.1, from
  https://github.com/aminforouzandehmehr/Piezo1-Channel-Model
  (SHA-256 `91406bb88352c18de4e5a7ec902411bf25ee7d639ea0e5ba3557ca2fcd95c208`).
  MIT License, Copyright (c) 2024 Amin Forouzandehmehr (`upstream/LICENSE_Piezo1`).
- Paper: Forouzandehmehr M, Ghosi S, Paci M, Hyttinen J, Koivumäki J. Mechanosensitive
  Channel Piezo1 in R403Q Hypertrophic Cardiomyopathy: A Computational Study.
  Computing in Cardiology 2023;50. DOI 10.22489/CinC.2023.359.
- The modified Piezo1 model of this release replaces the voltage dependent gate of
  version 1.1 by an availability (inactivation) gate with three parameters fitted to
  published channel recordings, adds a conductance scale, an adjustable reference
  stress, a first order filter of the native stress and a simplified Na+/Ca2+ ion
  accounting. The equations of version 1.1 remain available as an option
  (`X.gating = 'released'`, `X.accounting = 'published'`), checked by
  `verify_release.m` (check 1b).

## ChR2 (Williams et al. 2013)

- The ChR2 block of `model/hiMCES2026_ChR2_Piezo1.m` implements the four state
  ChR2(H134R) model published in: Williams JC, Xu J, Lu Z, Klimas A, Chen X, Ambrosi CM,
  Cohen IS, Entcheva E. Computational optogenetics: empirically-derived voltage- and
  light-sensitive channelrhodopsin-2 model. PLoS Comput Biol 2013;9(9):e1003220.
  DOI 10.1371/journal.pcbi.1003220.
- Equations and parameter values follow that paper. Where the paper rounds a
  constant (the prefactor of the recovery rate Gr and the photon flux constant), the
  value of the authors' MATLAB code (ModelDB accession 151549) is used; this is
  noted in the model file. No file or code of the authors is included in this release.

## Calibration data

The fitting targets and checks in `calibration/Piezo1_targets.m` are values read from
figures and text of the cited publications (Wu et al. 2017, Zheng et al. 2019,
Lewis et al. 2017, Albuisson et al. 2013), with the sources given at each value.

## Simulation results

The initial conditions and reference outputs were produced with the code of this
release (see `CHANGELOG.md`). They are covered by the MIT License in `LICENSE`.
