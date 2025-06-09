# ClinicalASL

ClinicalASL is an open-source processing tool for **Arterial Spin Labeling (ASL) MRI data** — notably for **multi-delay ASL and CVR mapping**, designed to support clinical research and reproducible imaging pipelines.

ClinicalASL generates multi-delay ASL derived images of **CBF, arrival time, and CVR**, ready for **PACS integration and clinical use**.  
It reads in **Philips multiframe DICOM** and exports derived images as **NIFTI** and **DICOM** (compatible with Philips multiframe).

It provides tools for data conversion, motion correction, parameter extraction, ASL quantification, visualization, and more.  
The project is currently transitioning from an original **MATLAB-based implementation** to a new, actively developed **Python version**.

---

## Project Maintainer

**Author:** Jeroen Siero  
**Institution:** University Medical Center Utrecht (UMCU), The Netherlands  
**Contact:** j.c.w.siero@umcutrecht.nl  
**Repository:** [https://github.com/JSIERO/ClinicalASL](https://github.com/JSIERO/ClinicalASL)

---

## Project Versions

| Version | Folder | Status | Notes |
|---------|--------|--------|-------|
| Python  | [python/](./python) | Actively developed | Primary version — supports modern pipelines and reproducibility |
| MATLAB  | [matlab/](./matlab) | Legacy (deprecated) | Provided for legacy users — no new features planned |

---

## Migration Notice

We are transitioning ClinicalASL to a fully Python-based pipeline to enable:

- Easier installation and portability
- Better integration with modern neuroimaging tools
- Improved maintainability and extensibility

**New users and contributors should start with the [python/](./python) version.**

The [matlab/](./matlab) version is preserved for compatibility with existing pipelines, but will not receive new features.

---

## Project Structure



```
ClinicalASL/
  python/      - Actively developed Python version
    run_pipeline.py           - Main entry point
    clinical_asl_pipeline/    - Python pipeline modules
    colormaps/                - Custom colormap .mat files (static resources)
    elastixfiles/             - Elastix parameter .txt files (static resources)
    requirements.txt          - Python dependencies
    README.md                 - Python-specific instructions
  matlab/      - Legacy MATLAB version (deprecated)
    *.m        - MATLAB scripts
    README.md  - MATLAB-specific instructions
  LICENSE      - BSD 3-Clause License
  README.md    - Project overview and instructions
```

## Python Version
3.12
### Installation

It is recommended to install ClinicalASL in a **Python virtual environment** (venv):

```bash
# Create venv (one-time)
cd ~/GITHUB/ClinicalASL
/usr/bin/python3 -m venv venv

# Activate venv
source venv/bin/activate

# Install required packages
pip install -r requirements.txt
```

### Usage

```bash
# From the project root:
source venv/bin/activate
python python/run_pipeline.py --input <data_dir> --output <results_dir> [other options]
```

### Dependencies

- numpy
- scipy
- matplotlib
- nibabel
- pydicom
- antspyx 

See [python/requirements.txt](./python/requirements.txt) for full list.

## MATLAB Version (Legacy)

The original MATLAB version of ClinicalASL is preserved in the [matlab/](./matlab) folder.

It is no longer actively developed and is provided only for compatibility with existing MATLAB-based pipelines.

## License

This project is licensed under the terms of the [BSD 3-Clause License](./LICENSE).

## Contributing

Contributions are welcome — especially to the Python version.

Please open an Issue or Discussion if you would like to contribute.

## Acknowledgments
ClinicalASL is developed and maintained by:

Jeroen Siero
University Medical Center Utrecht (UMCU), The Netherlands
j.c.w.siero@umcutrecht.nl.
Many thanks to Thomas Kirk for his help  and advice
## How to cite

If you use ClinicalASL in your research, please cite this repository or any related publications (to be added).
