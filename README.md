# hiMCES2026 with modified Piezo1 and ChR2

Simulations of optical and mechanical pacing in hiPSC-derived cardiomyocytes.
Version 1.0.0. Software DOI: [10.5281/zenodo.23430167](https://doi.org/10.5281/zenodo.23430167).

Source code: [aminforouzandehmehr/hiMCES2026-OptoSono](https://github.com/aminforouzandehmehr/hiMCES2026-OptoSono).

This is the model code of Fig. 8A and Supplementary Fig. S5 of:

Pettebone ME, Han JL, Forouzandehmehr A, Chesla O, Chua CJ, Pozo MR, Heinson YW, Liu W,
Entcheva E. Optogenetic and sonogenetic pacing of unmodified cardiac tissue upon quick
integration with specialized biological actuators. bioRxiv 2026.
DOI [10.64898/2026.10.19.760310](https://doi.org/10.64898/2026.10.19.760310).

The model is the published hiPSC-CM model hiMCES2026 (host equations and parameters
unchanged) with two added currents: ChR2 after Williams et al. 2013, and a modified
form of the Piezo1 model of Forouzandehmehr et al., Computing in Cardiology 2023.
Credits and licences of these components are in `NOTICE.md`.

## Requirements

MATLAB, base product only (no toolboxes). Tested with MATLAB R2026a Update 2 on macOS.

## Quick start

From this folder, in MATLAB:

```matlab
MasterCompute_OptoSono          % Fig. 8A and both controls (about 10 s)
FigureS5_Piezo1                 % Supplementary Fig. S5 (about 10 s)
ok = verify_release()           % equation checks and reproduction (about 20 s)
ok = verify_release('full')     % also re-conditions both cells and runs all channel checks (about 3 min)
```

Set `saveFigs = 1` in the two figure scripts to write PDF and PNG files to `results/`.
All scripts write only to `results/`; the supplied initial conditions and references are
never overwritten.

## Protocol of Fig. 8A

| Item | Value |
|---|---|
| Temperature, bath | 304 K (31 C), Tyrode solution (`classicVSoptic = 1`, `X.T_override = 304`) |
| Light (ChR2) | 10 ms pulses at 0.15 mW/mm2, at 0, 2 and 4 s (0.5 Hz) |
| Mechanical pulse (Piezo1) | 0.336 kPa effective stress at the Piezo1 gate, 0.5 s, from 5.0 s (1 s after the third light pulse) |
| Controls | the same cell without the mechanical pulse; the cell with the Piezo1 current off |
| Initial state | `initial_conditions/OptoSono_31C_Piezo1.mat` (40 states) and `..._Piezo1off.mat` (39 states), at the onset of a light pulse after 1000 s of 0.5 Hz light pacing |
| Solver | `ode15s`, RelTol 1e-7, AbsTol 1e-10, MaxStep 0.5 ms; integration cut at every pulse edge (`model/run_segments.m`) |

The mechanical input is an effective stress added at the Piezo1 activation gate. It is
not converted from, or calibrated to, an acoustic pressure or ultrasound intensity.
Only 0.5 s mechanical pulses were simulated; shorter pulses were not tested.

## Model

States 1 to 33: published hiMCES2026. States 34 and 35: Piezo1 activation gate and
availability. States 36 to 39: ChR2 (O1, O2, C2, p). State 40: native stress filter
(present only when `X.tau_trans` is set; not a Piezo1 gate). Fig. 8A uses all 40 states;
the Piezo1-off control has 39. State names and units: `initial_conditions/*.csv`.

Outputs `dati(1:48)` are those of the published model. Added outputs: 49 I_Piezo1 (A/F),
50 Ca2+ flux counted for Piezo1 (mM/s), 51 stress sensed by the activation gate (kPa),
52 steady state of the activation gate, 53 steady state of the version 1.1 voltage gate,
54 time constant of the second gate (s), 55 mechanical pulse on, 56 I_ChR2 (A/F),
57 light on, 58 ChR2 C1 fraction, 59 Na+ share of I_Piezo1 (A/F), 60 recovery time
constant (s).

Modified Piezo1 model (`X.gating = 'reduced'`):

    I_Piezo1 = g_P a h (Vm - E_P),   E_P = 0
    da/dt = (a_inf - a)/tau_a,       a_inf = 1/(1 + exp(-(5 sigma/sigma_ref - 2.5)/0.25))
    sigma = sigma_cell,s + sigma_stim,   d sigma_cell,s/dt = (sigma_cell - sigma_cell,s)/tau_tr
    dh/dt = (1 - h)/tau_rec - a h/tau_inact(v)
    tau_inact(v) = k_inact (0.3 ms) alpha(v) beta(v),  v = 1000 Vm
    alpha(v) = 1/sqrt(1 + exp((v - 60)/50)),  beta(v) = 900/(1 + exp((171 - v)/72))

| Parameter | Value in Fig. 8A | Basis |
|---|---|---|
| g_P | 91.17 S/F (300 x 0.3039 S/F; about 9 nS per cell) | conductance scale chosen for the simulation; not measured |
| sigma_ref | 0.336 kPa (half activation at 0.168 kPa) | chosen for the simulation; version 1.1 uses 0.056 kPa |
| tau_tr | 0.15 s | assumption |
| tau_a | 3.162 ms at 295 K | fitted (deactivation); activation assumed equal |
| k_inact | 0.9098 | fitted (voltage dependent decay) |
| tau_rec | 0.2023 s at 295 K | fitted (half recovery) |
| Q10 | 2 for tau_a, tau_inact, tau_rec from 295 K (factor 0.5359 at 304 K) | assumption |
| f_Ca | 0.05 of the inward current as Ca2+, rest Na+ | derived from permeability ratios; not measured |

`model/OptoSono_settings.m` returns these settings. The defaults of
`model/ChR2_Piezo1_opts.m` are the published values of version 1.1 and are not the
Fig. 8A configuration.

## Calibration and checks of the Piezo1 gates

The three kinetic parameters were fitted to the isolated channel (`model/Piezo1_channel.m`)
with the protocols of the experiments (`calibration/`):

- fitted: decay time constant at -60 to +20 mV (Wu et al. 2017; Zheng et al. 2019), 9 of
  10 within the plotted error bars (Zheng at +20 mV: 26.8 +/- 3.2 ms measured, 21.3 ms
  model); deactivation 3.86 ms (target 4.0 ms) and half recovery 0.144 s (target 0.140 s)
  (Lewis et al. 2017), within the assigned 20 % tolerance;
- checks not used in the fit, including the failed ones: human PIEZO1 decay at -80 mV
  8.33 ms (8.6 +/- 0.4 ms, Albuisson et al. 2013); current left after 300 ms 0.060
  (0.096 +/- 0.040); mouse decay at -80 mV 8.29 ms (11.9 +/- 0.6 ms, failed); decay at
  +40 to +100 mV 24.5 to 32.4 ms (33 to 95 ms, failed); slow phase of recovery and
  depression during trains of brief stimuli (failed).

`calibration/check_Piezo1.m` writes the full table (`reference_outputs/check_report.txt`
holds it for the frozen parameters). `calibration/fit_Piezo1.m` repeats the fit
(optional, several minutes); it writes to `results/` and does not change the frozen
parameters in `model/ChR2_Piezo1_opts.m`. The original fit report is in
`reference_outputs/fit_report.txt`.

## Limitations

- Single cell. The simulation does not include acoustic propagation, actuator spheroids
  or coupling to unmodified cardiac tissue, and does not give an acoustic threshold.
- Whole-cell myofilament stress is used as a proxy for membrane tension; its mapping to
  activation and the 0.15 s filter are not calibrated against tension measurements.
- The conductance scale, sigma_ref and tau_tr are settings or assumptions, not measurements.
- The two-gate form has one recovery time constant: it does not reproduce slow phases of
  recovery or the depression during trains of brief stimuli, and inactivates faster than
  measured above +20 mV.
- With E_P = 0 the Piezo1 current turns outward while Vm is positive during the
  mechanical pulse (up to +0.41 A/F in Fig. 8A); its size is not validated against
  cardiomyocyte recordings.
- Ion accounting is simplified (fixed 5 % Ca2+ share of inward Piezo1 current; outward
  current as K+ with fixed K_i). The ChR2 current enters the membrane potential but no
  ion balance. Calcium results are qualitative.
- Simulations are at 31 C; the experiments of the preprint were at room temperature.
- Only the variable-step double-precision solver above was used to produce and verify
  these results.

## Files

| Path | Content |
|---|---|
| `MasterCompute_OptoSono.m` | Fig. 8A and controls |
| `FigureS5_Piezo1.m` | Supplementary Fig. S5 and channel values |
| `condition_OptoSono.m` | conditioning of the two cells from `initial_conditions/start_31C_*.mat` |
| `verify_release.m` | equation checks and reproduction of the references |
| `model/` | `hiMCES2026_ChR2_Piezo1.m` (model), `ChR2_Piezo1_opts.m`, `OptoSono_settings.m`, `Piezo1_channel.m`, `run_segments.m`, `host_arguments.m`, `heav.m`, `passiveForces.m` |
| `analysis/` | upstroke times, APD90 and the Fig. 8A measurements |
| `calibration/` | targets with sources, channel protocols, checks, optional fit |
| `initial_conditions/` | conditioned states (MAT and CSV with state names and units) and conditioning start states |
| `reference_outputs/` | Fig. 8A traces (1 ms) and end states, measurements, threshold scan, Fig. S5 curve, channel measurements, check and fit reports |
| `upstream/` | unchanged published `hiMCES2026.m` and `Piezo1.m` (version 1.1) with their licences, used by `verify_release.m` |
| `verification/` | sampled host states for the equation checks |
| `MANIFEST.sha256` | SHA-256 checksums of all files |

## Licence and citation

MIT License (`LICENSE`); upstream notices in `NOTICE.md`. Please cite the software DOI
and the preprint (`CITATION.cff`). The ChR2 equations are those of Williams et al.,
PLoS Comput Biol 2013;9(9):e1003220 (DOI 10.1371/journal.pcbi.1003220); the host
model is Li et al., J Mol Cell Cardiol Plus 2026;15:100834
(DOI 10.1016/j.jmccpl.2026.100834); the original Piezo1 model is Forouzandehmehr et
al., Computing in Cardiology 2023;50 (DOI 10.22489/CinC.2023.359).

Zenodo DOI: [10.5281/zenodo.23430167](https://doi.org/10.5281/zenodo.23430167).
