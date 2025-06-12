# ClinicalASL Tool - Release v1.0.0

**Repository:** [https://github.com/JSIERO/ClinicalASL](https://github.com/JSIERO/ClinicalASL)  
**Author:** Jeroen Siero, UMC Utrecht  
**Contact:** j.c.w.siero@umcutrecht.nl  
**Release Date:** 2025-06-11  
**Version:** v1.0.0  

## Description

ClinicalASL is an automated processing tool for Arterial Spin Labeling (ASL) MRI for multidelay ASL, including:

- DICOM to NIfTI conversion
- Look-Locker correction
- Motion correction
- Quantitative perfusion parameter estimation:
    - Cerebral Blood Flow (CBF)
    - Arterial Arrival Time (AAT)
    - Arterial Transit Artefacts (ATA)
    - Cerebrovascular Reactivity (CVR)
- Export of results as NIfTI, PNG, and DICOM (PACS-ready) formats.

Supports baseline / stimulus ASL protocols (acetazolamide, hypercapnia, breath-hold, etc.)

## Intended Use

**This software is provided for post-processing of ASL MRI data for research and clinical evaluation purposes.**  
It is not CE marked or FDA cleared 

## Installation

See `INSTALL.md` or follow instructions in the GitHub repository.

## Configuration

Example config file: `config/config_default.json`  
You can supply your own `config.json` per acquisition protocol.

## Dependencies

- Python 3.11+
- HD-BET (brain extraction)
- ANTsPy or external ANTs
- See `requirements.txt`

## License

BSD 3-Clause License.

---

