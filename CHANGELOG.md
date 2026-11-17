
---

### `CHANGELOG.md`

```markdown
# CHANGELOG.md

ClinicalASL Pipeline — Changelog

---

## v1.1.0 — 2026-11-17

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