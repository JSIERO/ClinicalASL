# ClinicalASL Tool - Release v1.0.0

**Repository:** [https://github.com/JSIERO/ClinicalASL](https://github.com/JSIERO/ClinicalASL)  
**Author:** Jeroen Siero, UMC Utrecht  
**Contact:** j.c.w.siero@umcutrecht.nl  
**Release Date:** 2025-06-11  
**Version:** v1.0.0  

## Description

Official ClinicalASL processing pipeline for quantitative multi-delay ASL MRI data — designed for clinical use with PACS export and visual reporting.

- DICOM to NIfTI conversion
- Look-Locker correction
- Motion correction
- Quantitative perfusion parameter estimation:
    - Cerebral Blood Flow (CBF)
    - Arterial Arrival Time (AAT)
    - Arterial Transit Artefacts (ATA)
    - Cerebrovascular Reactivity (CVR)
- Export of results as NIfTI, PNG, and DICOM (PACS-ready) formats.

Supports baseline / stimulus ASL protocols (acetazolamide/diamox, hypercapnia, breath-hold, etc.)

## Intended Use

**This software is provided for post-processing of ASL MRI data for research and clinical evaluation purposes.**  
It is not CE marked or FDA cleared 

## Installation

See `INSTALL.md` or follow instructions in the GitHub repository.

## Configuration

Example config file containing scan and analysis parameters: `config/config_default.json`  
You can supply your own `config.json` per acquisition protocol.

example `config/config_default.json`:
```json
{
    "version": "v1.0-MRI_DIAMOX_MDLL_preACZ_postACZ-2025",
    "tau": 2,
    "N_BS": 4,
    "readout": "2D",
    "labeleff": 0.85,
    "lambda": 0.9,
    "T1t": 1.3,
    "T1b": 1.65,
    "FWHM": 6,
    "range_cbf": [0, 100],
    "range_cvr": [-50, 50],
    "range_AAT": [0, 3.5],
    "range_ATA": [0, 125],
    "inference_method": "basil",
    "ASL_CONTEXT": ["baseline", "stimulus"],
    "context_study_tags": ["preACZ", "postACZ"]
}
```
## Dependencies

- Python 3.11+
- HD-BET (brain extraction)
- ANTsPy or external ANTs
- See `requirements.txt`

## Installation
Please refer to INSTALL.md for full installation instructions.

Summary:
```bash
git clone https://github.com/JSIERO/ClinicalASL.git
cd ClinicalASL
conda create -n clinicalasl python=3.11 -y
conda activate clinicalasl
pip install -r requirements.txt
pip install git+https://github.com/MIC-DKFZ/HD-BET.git
```
Also install:
dcm2niix → https://github.com/rordenlab/dcm2niix
ANTs → via Conda or system package

## Running the pipeline
Example command line:

```bash
python run_pipeline.py /path/to/DICOM_INPUT /path/to/OUTPUT_FOLDER --config config/config_default.json
```
The pipeline will:
Save all intermediate files to the output folder

Generate:
NIfTI images for ASL derived images: CBF, AAT, CVR, ATA
PNG images for review
DICOM images for PACS export
Log file: clinicalasl.log

## License

BSD 3-Clause License.

---

