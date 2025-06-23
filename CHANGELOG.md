
---

### `CHANGELOG.md`

```markdown
# CHANGELOG.md

ClinicalASL Pipeline — Changelog

---

## v1.0.0 — June 2025

Initial official release of ClinicalASL:

- Modular ASL processing pipeline
- Baseline + stimulus ASL supported
- DICOM to NIfTI conversion with dcm2niix
- Look-Locker correction
- Motion correction with ANTs
- Quantification with sssv or vaby (BASIL)
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