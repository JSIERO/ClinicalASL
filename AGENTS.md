# ClinicalASL - notes for AI coding agents

Clinical multi-delay ASL pipeline (MRI DIAMOX, UMC Utrecht): baseline + stimulus (pre/post acetazolamide) ASL
DICOM in, CBF / AAT / ATA / CVR maps out as NIfTI, PNG and PACS-ready DICOM. Runs as one Docker image on the
IMAGR platform (`jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager`). Not CE marked / FDA cleared.
Sister project: FetalSVR (same container pattern, `..\FetalSVR`).

Detailed process files (pipeline internals, quantification, DICOM export, running/debugging, run comparison,
release): `.claude/agents/*.md` and `.claude/skills/*/SKILL.md`. Read the one that matches the task, not all.
Current state (open PRs, pending runs, release steps): `HANDOFF.md` if present (not in git). Read it at the start
of a session; overwrite it at the end when the state changed (~20 lines, no patient data).

## Layout
- `python/run_pipeline.py`: entry (Look-Locker multi-delay, `config_default.json`). Args: `inputdir outputdir workingdir
  [--mode cvr|single] [--inference-method ssvb|vaby] [--config cfg.json]`. `--mode single` (or config `mode`): one
  Look-Locker scan without ACZ -> CBF/AAT/ATA only (`ASL_CONTEXT` = baseline, no study-tag filter = newest SOURCE
  ASL series, no step 11, no CVR; set in `main_pipeline.apply_analysis_mode`). `--inference-method` overrides the
  config's `inference_method` (default `ssvb`) only when given.
- `python/clinical_asl_pipeline/main_pipeline.py` `mri_diamox_umcu_clinicalasl_cvr()`: the steps, per context
  (`ASL_CONTEXT` baseline/stimulus <-> `context_study_tags` preACZ/postACZ):
  1 paths, 2 dcm2niix (rename + convert; `asl_convert_dicom_to_nifti.py`), 3 newest SOURCE series per tag
  (`dicomseries_description_patterns`), 4 DICOM params, 5 Look-Locker correction, 6 control/label interleave,
  7 HD-BET brain mask (`utils/run_bet_mask.py`, CLI `hd-bet`), 8 ANTsPy motion correction (`random_seed=1`), 9 outlier removal,
  10 three QASL fits (`asl_qasl_analysis.py`, CLI `qasl`): all PLDs -> AAT, PLD 2..last -> CBF, PLD 1-2 artoff -> ATA,
  11 ANTsPy rigid stimulus->baseline (`random_seed=1`), 12 results + CVR (`asl_save_results_cbfaatcvr.py`).
- DICOM export: `utils/save_data_dicom_grayscale.py` / `save_data_dicom_color.py` (Enhanced MR multi-frame /
  MR Image Storage; RealWorldValueMapping from the map's own scaling, CVR stored as int16) on shared helpers in
  `utils/dicom_common.py`; `utils/save_png_to_dicom.py` (Secondary Capture, series 999, one SeriesDescription, map name
  in ImageComments). Colour series = grayscale series number + 100 (CVR 888 / 988). Any export failure fails the run.
- External tools via `utils/run_command_with_logging.py` (shell=True; on non-zero exit logs the last 40 stderr lines, raises).
- Logging: `utils/setup_logging.py` -> stdout + `<outputdir>/clinicalasl.log` (one line per DICOM series, no patient
  identifiers; fatal errors with traceback); `config_used.json` in outputdir.
- Version: `python/clinical_asl_pipeline/__version__.py` (single source; the Docker build overwrites it with APP_VERSION).
  Keep `"version"` in `config/config_default.json`, the README, the `load_parameters.py` fallback and the Dockerfile
  `APP_VERSION` default equal to it.
- `tools/`: run comparison and DICOM validation scripts (compare_maps, compare_all, compare_dicom, run_dciodvfy,
  make_preacz_input, nightly_verify.ps1). Host python has numpy/nibabel/pydicom.
- `matlab/`: original MATLAB implementation (+ Windows exe's). Reference only, not in the image; don't search it unless asked.
- `requirements.txt`: ClinicalASL deps on top of QASL. `antspyx==0.5.4` pinned (0.6.x registration is bad).

## Dependencies (all in the image, env `qi`, python 3.11)
- QASL 1.4.2 (Quantified Imaging, licensed, private bitbucket via SSH key, commits in `docker/qasl-commits.txt`).
  Key lives outside the repo: `G:\DATA\MOYAMOYA\IMAGER\qi-key-jeroen`. Never copy it into the repo, a layer, a log or chat.
- HD-BET v2.0.1 (CPU torch), parameters in `/opt/hd-bet_params`.
- dcm2niix v1.0.20220720 built with OpenJPEG + CharLS (PACS single-frame DICOM is JPEG 2000; conda-forge can't read it).
- Dockerfile details and the release checklist: `.claude/skills/docker-release/SKILL.md`. The old `:latest` image
  (26 GB) has the QASL key in a layer: never push it anywhere public.

## Build / run
```
docker build --secret id=qasl_key,src=G:\DATA\MOYAMOYA\IMAGER\qi-key-jeroen --build-arg APP_VERSION=1.2.1 -t jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:1.2.1 .
docker run --rm --network none -v G:\DATA\MOYAMOYA\IMAGER\ASLtest:/input:ro -v G:\DATA\MOYAMOYA\IMAGER\outputfolder:/output -v G:\DATA\MOYAMOYA\IMAGER\workingdir:/workingdir jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:1.2.1 /input /output /workingdir --inference-method ssvb
```
Working-tree code without rebuild: `docker run ... -w /fix -e PYTHONPATH=/fix -v <repo>\python:/fix:ro --entrypoint python
<image> /fix/run_pipeline.py ...` (or `/fix/export_results.py <workingdir> <outputdir>` = step 12 only, ~10 min).
Use PowerShell, not Git Bash, for docker runs with Windows paths.

## Test data and verification
- `G:\DATA\MOYAMOYA\IMAGER\ASLtest`: one PACS export (single-frame DICOM, NOT anonymized). Never in the repo or on
  ProtonDrive; never print patient tags; never read a whole log, grep it.
- Test runs: `G:\DATA\MOYAMOYA\IMAGER\tests\<name>\`; no new top-level folders in IMAGER. `IMAGER\dciodvfy\` holds the
  DICOM validation logs (final tool-validation step).
- Verification runs (full pipeline, export-only, dciodvfy, comparisons) are scheduled at night with
  `tools/nightly_verify.ps1` unless the user says otherwise; report the results the next session.
- No automated tests in the repo.

## Rules
- Never commit DICOM, NIfTI, logs or keys (see .gitignore). Never print patient tags (name, ID, birth date).
- IMAGR: logs to stdout/stderr, runs with `--network none`, non-zero exit on failure.
- Image tag = app version (`__version__`, set via APP_VERSION build-arg).
- Every change gets a new version X.Y.Z (QMS section 4 classes): minor (output-neutral: logging, docs, packaging,
  refactor) -> Z+1; major (changes quantitative output) -> Y+1, Z=0; very major (acquisition: sequence, protocol,
  scanner/coil) -> X+1, Y=Z=0. Bump all sync points above plus the CHANGELOG heading, and set the ClinicalASL QMS
  Procedure doc (claude.ai/artifact/JgBQkKFo6UGFfSCaVbRQ79) to the same version with a document-history row.
- Python only, concise code. Commit author JSIERO <76999501+JSIERO@users.noreply.github.com>.
- Log every user-visible change in `CHANGELOG.md` (newest first, short bullets).
