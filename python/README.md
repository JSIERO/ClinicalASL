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
    "ASL scan": "multi-delay Look-Locker",
    "tau": 2,
    "N_BS": 4,
    "readout": "2D",
    "labeleff": 0.85,
    "lambda": 0.9,
    "T1t": 1.3,
    "T1b": 1.65,
    "FWHM": 6,
    "outlier_factor": 2.5,
    "range_cbf": [0, 100],
    "range_cvr": [-50, 50],
    "range_AAT": [0, 3],
    "range_ATA": [0, 125],
    "inference_method": "vaby",
    "device": "cpu",
    "ASL_CONTEXT": ["baseline", "stimulus"],
    "context_study_tags": ["preACZ", "postACZ"]
}
```
## Dependencies

- Python 3.11+
- [QASL](https://quantified-imaging.com/) by Quantified Imaging (license key required)
- [dcm2niix](https://github.com/rordenlab/dcm2niix) (v1.0.20230411 or newer)
- [HD-BET](https://github.com/MIC-DKFZ/HD-BET) (for brain extraction)
- [ANTsPy](https://github.com/ANTsX/ANTsPy) (for image registration)
- Additional Python packages (see `requirements.txt`)

## Installation

Complete installation instructions are available in [INSTALL.md](INSTALL.md).

### Quick Start (QI Conda Environment)
```bash
# Clone repositories
git clone https://github.com/JSIERO/ClinicalASL.git
git clone https://bitbucket.org/quantified-imaging/qasl_setup.git

# Install QASL (requires license key)
cd qasl_setup
./qasl_setup --yes --key=<your-qasl-key>

# Set up environment
conda activate qi
conda install -n qi -c conda-forge pip dcm2niix

# Install HD-BET
git clone https://github.com/MIC-DKFZ/HD-BET.git
cd HD-BET
pip install -e .

# Install ClinicalASL requirements
cd ../ClinicalASL/python
pip install -r requirements.txt

## Running the pipeline
Example command line:

```bash
python run_pipeline.py /path/to/DICOM_INPUT /path/to/OUTPUT_FOLDER \
    --inference-method [vaby|ssvb] \
    --config /path/to/config.json
```

##Output Structure
The pipeline generates:

- NIfTI images: CBF, AAT, CVR, ATA maps
- Quality control PNGs: Visualizations for review
- DICOM files: PACS-compatible outputs
- Log file: clinicalasl.log with processing details
- Configuration copy: config_used.json

## License

BSD 3-Clause License. See LICENSE for details.
---

