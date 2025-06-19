# INSTALL.md

## ClinicalASL Installation Guide

This guide explains how to install the ClinicalASL processing pipeline.

---

## Table of Contents
1. [System Requirements](#system-requirements)
2. [Conda Environment Setup](#conda-environment-setup)
3. [QASL Installation](#qasl-installation)
4. [HD-BET Installation](#hd-bet-installation)  
5. [ANTsPy Installation](#antspy-installation)
6. [ClinicalASL Installation](#clinicalasl-installation)
7. [Docker Installation (Optional)](#docker-installation-optional)
8. [Verification](#verification)
9. [Troubleshooting](#troubleshooting)

## System Requirements
- **OS**: Linux (Ubuntu 22.04 recommended) or macOS
- **Memory**: Minimum 16GB RAM (32GB recommended for large datasets)
- **Storage**: 20GB free disk space
- **GPU**: Not required but recommended for HD-BET

## Conda Environment Setup

```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
source $HOME/miniconda/bin/activate

# Create environment
conda create -n qi python=3.11 -y
conda activate qi
```
## QASL Installation
```bash
# Get QASL (requires license key from https://quantified-imaging.com/)
git clone https://bitbucket.org/quantified-imaging/qasl_setup.git
cd qasl_setup
./qasl_setup --yes --key=YOUR_LICENSE_KEY_HERE
cd ..
```
## HD-BET Installation
```bash
git clone https://github.com/MIC-DKFZ/HD-BET.git
cd HD-BET
git checkout v2.0.1 
pip install -e .
cd ..
```
## ClinicalASL Installation
```bash
git clone https://github.com/JSIERO/ClinicalASL.git
cd ClinicalASL/python
pip install -r requirements.txt
```
### dcm2niix Installation
```bash
conda install -n qi -c conda-forge dcm2niix=1.0.20220720 -y
```
## Verification
```bash
python -c "import HD_BET; import ants; from clinical_asl_pipeline import main_pipeline; print('Installation successful')"
```
## Troubleshooting
### QASL Key Errors
- Ensure your license key is valid
- Contact support@quantified-imaging.com if issues persist
### HD-BET Model Download Issues
```bash
mkdir -p ~/.hd-bet
wget https://zenodo.org/record/14445620/files/release_v1.5.0.zip -P ~/.hd-bet
unzip ~/.hd-bet/release_v1.5.0.zip -d ~/.hd-bet
```
### ANTsPy Installation Problems
```bash
conda install -n qi -c conda-forge gcc=12.1.0 -y
conda install -n qi -c conda-forge antspy -y
```
For additional support, please open an issue on our GitHub repository.
