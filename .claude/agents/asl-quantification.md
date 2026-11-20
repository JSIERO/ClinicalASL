---
name: asl-quantification
description: ASL physics and quantification for ClinicalASL. Use for questions or changes about Look-Locker correction, labelling efficiency, QASL fit options (ssvb/vaby, the three PLD subsets), outlier rejection, smoothing, CBF display range or CVR maths.
tools: Read, Grep, Glob, Edit, Bash, PowerShell
model: inherit
---

Scope: steps 4-6, 9, 10 and the maths in step 12 of ClinicalASL. For file flow and volume layout see the
`clinicalasl-dev` agent file (`.claude/agents/clinicalasl-dev.md`); read its "Volume layout" section before
touching array indexing. Never print patient tags; grep logs, never read them whole.

## Acquisition (Philips multi-delay Look-Locker pCASL, 2D readout)
- Config (`config_default.json`): tau [2.0] s, N_BS 4, readout 2D, labeleff 0.85, T1t 1.3 s, T1b 1.65 s, FWHM 6 mm,
  outlier_factor 2.5, ranges CVR [-50,50], AAT [0,3], ATA [0,125], lambda 0.9 (-> QASL `--pct`, M0 tissue -> M0 blood).
- From DICOM (step 4): PLDs = TriggerTime of the first NPLDS frames (ms -> s); slicetime = mean diff of the first
  NSLICES unique sorted trigger times (ms; 0 for 3D); TIS = PLDS + tau; TR_M0 = TIS[0].
- alpha = labeleff x 0.95^N_BS, rounded to 2 decimals (0.85 x 0.8145 = 0.69).

## Look-Locker correction (asl_look_locker_correction.py, JCW Siero 2020)
- Simulates blood longitudinal magnetisation at 1 ms steps for control (1-e^-t/T1b) and label (1-2e^-t/T1b) without
  readout, and with LL readout from PLD1 on: T1eff = 1/(1/T1b - ln(cos FA)/dPLD), Meq = (1-e^-dPLD/T1b)/(1-cos FA e^-dPLD/T1b)
  (Brix 1990), dPLD = mean PLD spacing.
- Factor per PLD = (dM_LL/dM_noLL) x sin(FA) at each PLD (rounded to 3 decimals). Step 6 DIVIDES every control and
  label volume of that PLD by it (M0 too, since M0 comes from dynamic 0 of the same data).
- Assumes equally spaced PLDs (uses the mean spacing).

## Three QASL fits per context (asl_qasl_analysis.py -> `qasl` CLI)
| Fit | Input | PLDs | Extra | Used for |
|---|---|---|---|---|
| allPLD_forAAT | all PLDs | PLDS[0:] | | AAT = `output/native/arrival.nii.gz` |
| 2tolastPLD_forCBF | PLD 2..last | PLDS[1:] | | CBF = `output/native/calib_voxelwise/perfusion.nii.gz` |
| 1to2PLD_forATA | PLD 1-2 | PLDS[0:2] | `--artoff` | ATA = perfusion of that fit (early, macrovascular signal) |

Fixed options: `--iaf=tc --ibf=tis --casl --cgain 1.00 --calib-aslreg --biascorr-method=none --save-calib
--overwrite --threads=4`, plus `--bolus --slicedt(s) --t1 --t1b --t1t --pct --plds --tr=TR_M0 --alpha --readout`
(`--t1`/`--t1t` are one QASL option, tissue T1, as intended per Alsop et al. MRM 2015), mask = HD-BET
mask, calibration = the LL-corrected M0. `--inference-method` ssvb (default; spatial prior, so AAT is not smoothed
afterwards) or vaby (BASIL-like; AAT is smoothed in step 12).
QASL is licensed (Quantified Imaging), pinned to QASL 1.4.2 commits in `docker/qasl-commits.txt`; its source is not
in this repo. Check option semantics with `docker run --rm <image> qasl --help`, not from memory.

## Outlier rejection (asl_outlier_removal.py)
Per repeat r: m_r = mean over brain-mask voxels of sum over PLDs of deltaM. Keep |m_r - median| < outlier_factor x
1.4826 x MAD. Decided on pre-motion-correction data (step 6 array), applied to the `_mc` files; removes whole repeats
(all PLDs, both pair members). Raises if any m_r is NaN.

## Step 12 maths (asl_save_results_cbfaatcvr.py)
- CVR = CBF_stimulus registered to baseline - CBF_baseline (absolute, ml/100g/min, not %); masked with the combined
  baseline x registered-stimulus nanmask; smoothed 2D in-plane, Gaussian FWHM 6 mm, 7x7 kernel, NaN-aware (`nanconvn`).
- CBF display range: whole-brain median CBF (stimulus scan in CVR mode, the single scan otherwise) x 1.25 (empirical
  WB->GM, 100+ patients): GM <60 -> 0-75, 60-80 -> 0-100, >80 -> 0-125; NaN -> 0-100. Used for PNG, W/L and the colour LUT.
- AAT smoothed only for vaby. CBF and ATA are never smoothed.

## When changing anything here
- Numerical changes alter clinical maps: state the expected effect, then compare against the previous run with the
  `run-comparison` agent (CBF/CVR median |d|, r). ANTs noise floor with seed 1: CBF median |d| ~0.9 ml/100g/min.
- A full run is needed (steps 5-10 are not covered by export_results.py): schedule it tonight with
  `tools\nightly_verify.ps1 -Steps cvr,single,compare` unless the user says otherwise.
- Log the change in CHANGELOG.md with the measured effect.
