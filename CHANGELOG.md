
---

### `CHANGELOG.md`

```markdown
# CHANGELOG.md

ClinicalASL Pipeline — Changelog

---

## v1.2.1 — 2026-11-19 (minor; verification in the nightly run, image not built yet)

- Logging: one line per DICOM series instead of five header lines per file; no PatientID, StudyID, StudyDate/Time
  or StudyInstanceUID in `clinicalasl.log` any more. Step 2 logs one line per skipped SeriesDescription instead of
  every skipped file path. A CVR log goes from ~1700 to ~300 lines.
- Failures: a fatal error is logged with its traceback in `clinicalasl.log`; a failing external command (dcm2niix,
  hd-bet, qasl) logs its last 40 stderr lines.
- DICOM/PNG export errors now stop the run with a non-zero exit (were logged and skipped, exit 0, series missing).
  An all-zero map is stored as zeros with a warning instead of failing on the scaling.
- QASL: config `lambda` is now passed as `--pct` (partition coefficient). Same value as the QASL default (0.9), so
  results do not change.
- Code: shared DICOM export code in `utils/dicom_common.py` (grayscale/colour modules 956 -> 507 lines incl. the
  shared module), same output intended; removed unused `asl_t1_from_m0.py`; `load_parameters` default config path fixed.
- Repo: `tools/` with run comparison and DICOM validation scripts (`compare_dicom.py`, `nightly_verify.ps1`, ...;
  `compare_maps.py` also reports median |d|, the QMS V4 criterion);
  `AGENTS.md` (instructions for AI coding tools, CLAUDE.md imports it); `.gitignore` excludes DICOM/NIfTI.

---

## v1.2.0 — 2026-11-17

- New `--mode single` (Docker/CLI, also `export_results.py`; config key `mode`): processes one multi-delay
  Look-Locker scan acquired without acetazolamide -> CBF, AAT, ATA (NIfTI, PNG, DICOM), no registration, no CVR.
  Uses the newest SOURCE ASL series; series/PNG names without preACZ/postACZ (e.g. "ASL CBF"). Default `--mode cvr`
  is unchanged.
  Image validated on ASLtest: single mode on the pre-ACZ scan (28.6 min, 99 DICOMs in 7 series) and CVR mode on the
  full exam (63 min, same 231 DICOMs in 15 series as 1.1.0; masks identical, map differences within ANTs
  run-to-run noise, e.g. CBF median |d| 0.9 ml/100g/min). dciodvfy errors only PatientOrientation on colour series.
- Config: new key `"mode": "cvr"`; config `"version"` is now the app version ("1.2.0", was
  "v1.1-MRI_DIAMOX_MDLL_preACZ_postACZ-2025").

---

## v1.1.0 — 2026-11-17

- ANTs motion correction and the stimulus->baseline rigid registration use a fixed random seed (`random_seed=1`).
  On ASLtest, repeat runs differ much less: motion-corrected data mean |d| 0.35 instead of 0.90 (signal ~334),
  registration 0.006 instead of 3.2. Not fully identical (multithreaded ITK); 1 thread would be, but ~3x slower.
- Brain mask: HD-BET on M0 + all control/label (see Fix below). On the test case the mask is ~1.5% larger (a
  1-voxel cortical rim, Dice 0.99); CBF/CVR agree with the M0-only mask at r 0.995/0.994.
- DICOM (dciodvfy clean of errors except PatientOrientation on the colour series): RealWorldValueMapping now has
  LUTLabel; colour series ConversionType = WSD; PNG series 999 has Laterality (as the map series), one
  SeriesDescription "ASL results (ClinicalASL, NOT CE MARKED)" for the whole series, map name in ImageComments.
- Colour CVR DICOM series number 988 (grayscale CVR 888 + 100, like the other colour series; was 888 for both).
- No author name in any SeriesDescription (PNG series 999 had "(ClinicalASL-Siero)"); header fields and log unchanged.
- New `python/export_results.py`: re-runs only the results/export step from an existing working directory, to
  test PNG/DICOM changes without a full run.
- Docker: lean multi-stage image (6.6 GB, CPU-only torch, QASL 1.4.2 pinned to the production commits).
- Removed the variable-TR (vTR) variant: `run_pipeline_vTR.py`, `main_pipeline_vTR.py`, `config_vTR.json`,
  the vTR branches in parameter extraction, data preparation and QASL, and the vTR MATLAB scripts/exe.
  Only multi-delay Look-Locker ASL is supported; other `ASL scan` types now raise an error.
- Outputs now state "ClinicalASL, NOT CE MARKED" (tool name + disclaimer; replaces "RESEARCH ONLY, NOT CE MARKED"):
  PNG titles, DICOM SeriesDescription (longest 50 of 64 characters), ImageComments and DerivationDescription.
- Fix: HD-BET brain mask now runs on the combined M0 + all control/label image as intended (it ran on M0
  alone; the combined image was written and deleted unused). Changes the brain mask and so the masked maps.
- Fix: `--inference-method` no longer defaults to `ssvb` in argparse, so `inference_method` from the config is
  used unless the option is given.
- Docs: README/INSTALL now show the `workingdir` argument, dcm2niix v1.0.20220720 built with JPEG 2000 support
  and the Docker build/run.

---

## v1.0.0 — June 2025

Initial official release of ClinicalASL:

- Modular ASL processing pipeline
- Baseline + stimulus ASL supported
- DICOM to NIfTI conversion with dcm2niix
- Look-Locker correction
- Motion correction with ANTs
- Quantification with ssvb or vaby (BASIL)
- Registration of post-stimulus to baseline space
- Computation of CVR maps
- PNG and DICOM export for PACS systems
- Configurable via `config/config_default.json`
- Logging to `clinicalasl.log`

---

## Planned for next versions

- Improved GUI wrapper (ClinicalASL-Viewer)
- Support multi-band ASL sequences
- BIDS input/output compatibility
- Auto-parameter tuning
- Batch processing scripts

---

Notes

Release author:  
Jeroen Siero  
UMC Utrecht  
j.c.w.siero@umcutrecht.nl

Repository: https://github.com/JSIERO/ClinicalASL