# INSTALL.md

## ClinicalASL Installation Guide

This guide explains how to install the ClinicalASL processing pipeline.

---

## 1️⃣ Prerequisites

- Python 3.11 (recommended via Conda)
- Git
- ANTs (Advanced Normalization Tools) → must be in PATH
- HD-BET (will be installed by pip)
- Linux or Windows with WSL / Ubuntu preferred
- dcm2niix installation:
    - Recommended: conda install -c conda-forge dcm2niix
    - Advanced: build from source → https://github.com/rordenlab/dcm2niix
    
## Recommended environment (Conda)

```bash
# Create and activate environment
conda create -n clinicalasl python=3.11 -y
conda activate clinicalasl

# Clone the repository
git clone https://github.com/JSIERO/ClinicalASL.git
cd ClinicalASL

# Install required Python packages
pip install -r requirements.txt


## 2️⃣ Clone the repository

git clone https://github.com/JSIERO/ClinicalASL.git
cd ClinicalASL

