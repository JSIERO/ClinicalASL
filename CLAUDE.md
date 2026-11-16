# ClinicalASL - notes for Claude

Clinical multi-delay ASL pipeline (MRI DIAMOX, UMC Utrecht): baseline + stimulus (pre/post acetazolamide) ASL
DICOM in, CBF / AAT / ATA / CVR maps out as NIfTI, PNG and PACS-ready DICOM. Runs as one Docker image on the
IMAGR platform (`jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager`). Not CE marked / FDA cleared.
Sister project: FetalSVR (same container pattern, `..\FetalSVR`).

## Layout
- `python/run_pipeline.py`: entry (Look-Locker multi-delay, `config_default.json`). Args: `inputdir outputdir workingdir
  [--inference-method ssvb|vaby] [--config cfg.json]`. `--inference-method` defaults to `ssvb` in argparse, so the
  config value is never used unless the default is changed.
- `python/run_pipeline_vTR.py`: variable-TR 3D variant (`config_vTR.json`, `main_pipeline_vTR.py`). BROKEN as entry:
  `main()` calls `run_pipeline()` without `workingdir`. The QASL step in `main_pipeline_vTR.py` is commented out.
- `python/clinical_asl_pipeline/main_pipeline.py` `mri_diamox_umcu_clinicalasl_cvr()`: the steps, per context
  (`ASL_CONTEXT` baseline/stimulus <-> `context_study_tags` preACZ/postACZ):
  1 paths, 2 dcm2niix (rename + convert; `asl_convert_dicom_to_nifti.py`), 3 newest SOURCE series per tag
  (`dicomseries_description_patterns`), 4 DICOM params, 5 Look-Locker correction, 6 control/label interleave,
  7 HD-BET brain mask (`utils/run_bet_mask.py`, CLI `hd-bet`), 8 ANTsPy motion correction, 9 outlier removal,
  10 three QASL fits (`asl_qasl_analysis.py`, CLI `qasl`): all PLDs -> AAT, PLD 2..last -> CBF, PLD 1-2 artoff -> ATA,
  11 ANTsPy rigid stimulus->baseline, 12 results + CVR (`asl_save_results_cbfaatcvr.py`).
- DICOM export: `utils/save_data_dicom_grayscale.py` / `save_data_dicom_color.py` (Enhanced MR multi-frame /
  MR Image Storage; RealWorldValueMapping set from the derived map's own scaling, CVR stored as int16),
  `utils/save_png_to_dicom.py` (Secondary Capture, series 999). Colour series = grayscale series number + 100.
- External tools are called through `utils/run_command_with_logging.py` (shell=True, raises on non-zero exit).
- Logging: `utils/setup_logging.py` -> stdout + `<outputdir>/clinicalasl.log`; `config_used.json` in outputdir.
- Version: `python/clinical_asl_pipeline/__version__.py` (single source; the Docker build overwrites it with APP_VERSION).
- `matlab/`: original MATLAB implementation (+ Windows exe's elastix, bet2, dcm2niix). Reference only, not in the image.
- `requirements.txt`: ClinicalASL deps on top of QASL. `antspyx==0.5.4` pinned (0.6.x registration is bad).

## Dependencies (all in the image, env `qi`, python 3.11)
- QASL (Quantified Imaging, licensed): `qasl_setup` pip-installs qasl, ssvb, vaby, toblerone, regtricks, ... from
  private bitbucket repos via SSH key. ssvb/vaby use TensorFlow; toblerone pulls vtk. Key lives outside the repo:
  `G:\DATA\MOYAMOYA\IMAGER\qi-key-jeroen`. Never copy it into the repo, a layer, a log or chat.
- HD-BET v2.0.1 (torch, nnunetv2), parameters from zenodo `release_v1.5.0.zip` in `/opt/hd-bet_params`
  (`HD_BET/paths.py` patched).
- dcm2niix v1.0.20220720 built with OpenJPEG + CharLS: PACS/IMAGR single-frame DICOM is JPEG 2000 compressed;
  conda-forge dcm2niix cannot read it.

## Docker
- `Dockerfile` (lean, multi-stage like FetalSVR): stage `dcm2niix` builds dcm2niix; stage `pyenv` makes micromamba env
  `/opt/conda/envs/qi`, installs CPU-only torch FIRST (so no nvidia-*/triton CUDA wheels, ~3.4 GB), then runs
  qasl_setup with the key as BuildKit secret, HD-BET, requirements.txt, and strips test/dev packages, caches, headers.
  Final = ubuntu:22.04 + runtime libs (libgomp1, libopenjp2-7, libcharls2) + env + params + `python/`. A FINAL CHECK
  RUN fails the build if imports, CPU torch, qasl/hd-bet CLIs or dcm2niix are missing.
- The old image (`G:\DATA\MOYAMOYA\IMAGER\Dockerfile`, tag `latest`, 26 GB) cloned the repo from GitHub, kept
  compilers, CUDA torch, 3.7 GB pip cache, and COPY'd the QASL key into a layer (still in its history: do not push
  that image anywhere public).
- `docker-entrypoint.sh`: `bash|python|qasl|hd-bet|dcm2niix|...` pass through, otherwise `run_pipeline.py "$@"`.
  Kept LF via `.gitattributes` (core.autocrlf=true on this machine); the build also strips CR.

## Build / run
```
docker build --secret id=qasl_key,src=G:\DATA\MOYAMOYA\IMAGER\qi-key-jeroen --build-arg APP_VERSION=1.0.4 -t jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:1.0.4 .
docker run --rm --network none -v G:\DATA\MOYAMOYA\IMAGER\ASLtest:/input:ro -v G:\DATA\MOYAMOYA\IMAGER\outputfolder:/output -v G:\DATA\MOYAMOYA\IMAGER\workingdir:/workingdir jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:1.0.4 /input /output /workingdir --inference-method ssvb
```
Debug shell: `docker run --rm -it ... <image> bash`. On Linux hosts (RIBS) the output folder must be writable by the
container user (chmod 777). Use PowerShell, not Git Bash, for docker runs with Windows paths.

## Test data
- `G:\DATA\MOYAMOYA\IMAGER\ASLtest`: one PACS export (single-frame DICOM, not anonymized). Never in the repo or on
  ProtonDrive; never echo patient tags or names. Outputs to `outputfolder\`, intermediates to `workingdir\`.
- No automated tests in the repo.

## Known issues (open)
- `utils/run_bet_mask.py`: `cmd` is built before `inputdata_path` is switched to the combined M0+ASL image, so HD-BET
  always runs on M0 alone; the combined `temp_for_bet.nii.gz` is written and deleted unused.
- `run_pipeline_vTR.py` missing `workingdir` argument (see above).
- README/INSTALL mention `python run_pipeline.py in out` without workingdir and log `clinicalasl.log`; dcm2niix
  version in README (1.0.20230411+) differs from the one actually used (1.0.20220720).

## Rules
- Never commit DICOM, NIfTI, logs or keys (see .gitignore). Never print patient tags (name, ID, birth date).
- IMAGR: logs to stdout/stderr, runs with `--network none`, non-zero exit on failure.
- Image tag = app version (`__version__`, set via APP_VERSION build-arg).
- Python only, concise code. Commit author JSIERO <76999501+JSIERO@users.noreply.github.com>.
- Log every user-visible change in `CHANGELOG.md` (newest first, short bullets).
