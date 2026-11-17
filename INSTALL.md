# INSTALL.md

## ClinicalASL Installation Guide

This guide explains how to install the ClinicalASL processing pipeline.

---

## Table of Contents
1. [System Requirements](#system-requirements)
2. [Conda Environment Setup](#conda-environment-setup)
3. [dcm2niix installation](#dcm2niix-installation)
4. [ClinicalASL Installation](#clinicalasl-installation)
5. [QASL Installation](#qasl-installation)
6. [HD-BET Installation](#hd-bet-installation)  
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

# Create environment (important!)
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
## ClinicalASL Installation
```bash
git clone https://github.com/JSIERO/ClinicalASL.git
cd ClinicalASL
pip install -r requirements.txt
```

## HD-BET Installation
```bash
git clone https://github.com/MIC-DKFZ/HD-BET.git
cd HD-BET
git checkout v2.0.1 
pip install -e .
cd ..
```

## dcm2niix Installation
Use v1.0.20220720 built with JPEG 2000 support (OpenJPEG + CharLS), as in the Dockerfile. PACS single-frame DICOM is JPEG 2000 compressed and the conda-forge build cannot read it.
```bash
sudo apt-get install -y git cmake make g++ pkg-config zlib1g-dev libopenjp2-7-dev libcharls-dev
git clone --depth 1 --branch v1.0.20220720 https://github.com/rordenlab/dcm2niix.git
cd dcm2niix && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OPENJPEG=ON -DUSE_CHARLS=ON ..
make -j && sudo cp bin/dcm2niix /usr/local/bin/
cd ../..
```

## Docker Installation (Optional)
The image contains QASL, HD-BET (with parameters), dcm2niix and ClinicalASL. The QASL key is passed as a BuildKit secret and is never stored in the image.
```bash
docker build --secret id=qasl_key,src=/path/to/qasl-key --build-arg APP_VERSION=1.2.0 -t jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:1.2.0 .
docker run --rm --network none -v /path/to/input:/input:ro -v /path/to/output:/output -v /path/to/workingdir:/workingdir jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:1.2.0 /input /output /workingdir
```
## Verification
```bash
python -c "import HD_BET; import ants; from clinical_asl_pipeline import main_pipeline; print('Installation successful')"
```
## Troubleshooting
### DCM2NIIX conversion of single-frame DICOMS from PACS
-   dcm2niix needs JPEG 2000 support: build v1.0.20220720 from source as described above

### QASL Key Errors
- Ensure your license key is valid
- Contact support@quantified-imaging.com if issues persist

### HD-BET Model Download Issues (looks default in ~/hd-bet folder)
```bash
mkdir -p ~/.hd-bet
wget https://zenodo.org/record/14445620/files/release_v1.5.0.zip -P ~/.hd-bet
unzip ~/.hd-bet/release_v1.5.0.zip -d ~/.hd-bet
```
For additional support, please open an issue on our GitHub repository.
