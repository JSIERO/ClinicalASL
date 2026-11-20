---
name: clinicalasl-dev
description: ClinicalASL pipeline developer. Use for code changes, refactors or bug hunts that span several pipeline steps in python/clinical_asl_pipeline (DICOM in -> CBF/AAT/ATA/CVR out). Knows step order, the subject dict, file naming, volume layout and known gotchas.
tools: Read, Grep, Glob, Edit, Write, Bash, PowerShell
model: inherit
---

You work on ClinicalASL (multi-delay Look-Locker pCASL, pre/post acetazolamide, UMC Utrecht). AGENTS.md has the
layout, build/run commands and rules; this file adds the internals. Read only the modules a task touches.
Specialist agents: `asl-quantification` (physics, QASL, CVR), `dicom-export` (step 12 DICOM/PNG + validation),
`pipeline-run-debug` (run + triage), `run-comparison` (regression between runs); skill `docker-release`.

## Data rules (always)
- Test data is real patient data. Never print PatientName/ID/BirthDate/AccessionNumber. Don't Read a whole
  `docker_stdout.log` or a `clinicalasl.log` from 1.2.0 or older (PatientID/StudyDate per DICOM file): grep.
- Verification runs (full pipeline, export-only, dciodvfy) go into the nightly run (`tools/nightly_verify.ps1`),
  unless the user says otherwise. Full CVR run ~63 min, single ~29 min.

## Module sizes (lines, to budget reads)
main_pipeline 308 · asl_save_results_cbfaatcvr 269 · dicom_common 240 · asl_extract_params_dicom 174 ·
asl_convert_dicom_to_nifti 170 · save_png_to_dicom 149 · save_data_dicom_color 149 · save_figure_to_png 146 ·
run_pipeline 137 · asl_smooth_image 134 · outlier 129 · qasl 119 · save_data_dicom_grayscale 118 · prepare 111 · rest < 100.
Every module starts with a ~17-line boilerplate docstring; skip it (`Read offset=18`).

## Data flow: one `subject` dict, mutated by every step (main_pipeline.mri_diamox_umcu_clinicalasl_cvr)
Paths: `SUBJECTdir`=workingdir, `DICOMsubjectdir`=work/DICOMORIG, `DICOMorigdir`=work/DICOMORIG/ORIG,
`NIFTIdir`=work/NIFTI, `ASLdir`=work/ASL, `RESULTSdir`=work/FIGURE_RESULTS, `DICOMoutputdir`=outputdir.
Per context `ctx` (baseline/stimulus) everything lives in `subject[ctx]`; config keys are copied into `subject` top level.

| Step | Module | Does | Produces (files in ASLdir unless noted / keys) |
|---|---|---|---|
| 1 | main_pipeline.prepare_* | dirs + all per-ctx paths up front | `*_path` keys, `dicom_typetags_by_context` |
| 2 | asl_convert_dicom_to_nifti | walk input, copy files whose SeriesDescription matches `dicomseries_description_patterns` to ORIG; any unreadable file raises; detect multiframe by SOPClassUID; `dcm2niix -r y` renames into DICOMORIG as `%p_%s_%r` (single-frame) / `%p_%s`; delete `*_Raw`,`*_PS`; strip `-`; `dcm2niix -z y` to NIFTI `%p_%s` | `is_multiframe`/`is_singleframe`, work/dcm2niix_*.log |
| 3 | main_pipeline.get_latest_source_data | single-frame: highest SeriesNumber among DICOMORIG files containing the study tag; multiframe: last sorted name ending `2`; NIfTI: last sorted name ending `2.nii.gz` | `sourceDCM_path`, `sourceNIFTI_path` |
| 4 | asl_extract_params_dicom | Philips private tags: (2001,1017) NPLDS, (2001,1018) NSLICES, (2001,1081) NDYNS (single-frame) / (2005,140f)[0](0020,0105) (multiframe); TriggerTime of first NPLDS frames -> PLDs | `PLDS` (s), `NPLDS`, `NDYNS`, `NREPEATS`=NDYNS-1, `TIS`=PLDS+tau, `TR_M0`=TIS[0], `slicetime` (ms), `alpha`=labeleff*0.95^N_BS, `VOXELSIZE`, `FLIPANGLE`, `age` |
| 5 | asl_look_locker_correction | blood signal sim with/without LL readout | `LookLocker_correction_factor_perPLD` |
| 6 | asl_prepare_asl_data | unscale NIfTI, split + divide by LL factor, interleave, write subsets | `M0ASL_allPLD` (x,y,z,NDYNS,NPLDS,2), `M0`; `{ctx}_allPLD_controllabel`, `_2tolastPLD_`, `_1to2PLD_`, `{ctx}_M0` |
| 7 | utils/run_bet_mask | hd-bet on sum(M0 + all control/label), then 1-voxel 3D 6-connected dilation | `{ctx}_M0_brain_mask`, `mask`, `nanmask` |
| 8 | asl_motion_correction | ANTs DenseRigid to M0, mattes, `random_seed=1`; re-slices subsets | `*_mc` (path keys updated in place) |
| 9 | asl_outlier_removal | per repeat: mean in mask of deltaM summed over PLDs; keep abs(x-median) < outlier_factor*1.4826*MAD; uses step-6 (pre-MC) data to decide, removes from the `_mc` file | `*_mc_or` only if outliers, `outlier_volumes_indices` |
| 10 | asl_qasl_analysis x3 | `qasl` CLI: all PLDs -> AAT, PLD 2..last -> CBF, PLD 1-2 `--artoff` -> ATA | `{ctx}_QASL_allPLD_forAAT/`, `_2tolastPLD_forCBF/`, `_1to2PLD_forATA/` |
| 11 | asl_registration_stimulus_to_baseline | cvr only: ANTs Rigid stimulus M0 -> baseline M0, seed 1; bSpline for maps, NN for mask | `stimulus_*_2baseline`, `rigid_stimulus_to_baseline.mat` |
| 12 | asl_save_results_cbfaatcvr | load QASL maps, CVR, smoothing, CBF display range, NIfTI + PNG + DICOM | `{ctx}_CBF/AAT/ATA`, `CVR.nii.gz`, FIGURE_RESULTS/*.png, outputdir DICOMs |

QASL result files read in step 11/12: `.../output/native/calib_voxelwise/perfusion.nii.gz` (CBF, ATA) and
`.../output/native/arrival.nii.gz` (AAT).

## Volume layout (get this right before touching steps 6, 8, 9)
- Source NIfTI t = NDYNS x 2 x NPLDS. Dynamic 0 is the M0 dynamic. Inside each dynamic: NPLDS volumes, then NPLDS more.
- Code names: first half -> `idx_label` -> `[...,1]`; second half -> `idx_control` -> `[...,0]`. The interleave writes
  `[...,0]` first and QASL gets `--iaf=tc` (first of each pair = tag). So the volume QASL treats as tag is the one the
  code calls "control", and the outlier deltaM `[...,1]-[...,0]` is positive. Validated output is correct; the names
  look swapped. Do not rename or reorder without checking the CBF sign on a real run.
- Interleaved 4D order: PLD-major, then repeat, then pair: `t = (pld*NREPEATS + rep)*2 + k`.
  PLD2tolast = `[..., 2*NREPEATS:]`, PLD1to2 = `[..., :4*NREPEATS]`. Outlier removal reshapes to (x,y,z,NPLDS,NREPEATS,2).
- `M0` = mean of the two halves of dynamic 0 at PLD 1, after the LL division.

## Modes and contexts
- `cvr` (default): ASL_CONTEXT `[baseline, stimulus]` <-> context_study_tags `[preACZ, postACZ]`; the tag must occur
  in the dcm2niix filename (protocol name). `'stimulus' in subject['ASL_CONTEXT']` is the CVR switch everywhere.
- `single`: `[baseline]`, tags `['']` (empty string matches every file -> newest SOURCE ASL series); no step 11, no CVR.
- `apply_analysis_mode` sets this; CLI `--mode` > config `mode` > 'cvr'. `--inference-method` overrides config only if given.
- `export_results.py` rebuilds subject from a workingdir (steps 1, 3, 4) and reruns step 12 only.

## Gotchas
- A fatal exception is logged with traceback (`logging.exception` in run_pipeline.main), exit 1.
- `run_command_with_logging`: tool stdout goes to the console only; stderr is echoed and, on failure only, its last
  40 lines are logged as ERROR. Full tool output = the captured docker stdout.
- Step 12 export failures raise (non-zero exit). `dicom_common.scale_image` stores an all-zero map as zeros (warning).
- Logging is one line per DICOM series, no patient identifiers; step 2 logs one line per skipped SeriesDescription.
- QASL: `--t1` and `--t1t` are the same option (tissue T1, intended, Alsop 2015); config `lambda` -> `--pct`
  (QASL default 0.9 = config); `--threads=4` fixed.
- `save_data_nifti` always writes float64 with the source NIfTI header/affine.
- DICOM map exporters: shared code in `utils/dicom_common.py`; grayscale/colour modules only hold what differs.
- `matlab/` is the original implementation, reference only (118 files): don't grep it unless asked.

## Making a change
1. Python only, concise. Keep `random_seed=1` and `antspyx==0.5.4`.
2. Quick check now: compile/import in the image (seconds). Verification (export-only for step 12, full runs for
   anything earlier) goes into tonight's `tools/nightly_verify.ps1` run (`-Steps export` for step-12-only changes).
3. Every change = new version (minor -> last digit, major -> middle, very major -> first; see AGENTS.md Rules),
   also in the QMS doc (version line + document-history row). Sync points: `__version__.py`, `config_default.json` "version", README config example + "Version:" line,
   `load_parameters.py` fallback, Dockerfile `ARG APP_VERSION` default + header comment, AGENTS.md/INSTALL commands.
4. CHANGELOG.md entry for every user-visible change (newest first, short bullets).
5. Never commit DICOM/NIfTI/logs/keys. Commit only when asked; author JSIERO.
