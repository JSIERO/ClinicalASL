#!/usr/bin/env python3
"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Main pipeline script for processing ASL MRI data.
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMC Utrecht (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    This script runs the main pipeline for processing ASL data, including DICOM to NIfTI conversion,
    parameter extraction, Look-Locker correction, data preparation, brain extraction, motion correction,
    quantification, registration, and saving of CBF, AAT, ATA, and CVR results for a subject.
    Results are saved as NIfTI, PNG, and DICOM files for further analysis and PACS export.

License: BSD 3-Clause License
"""
import argparse
from clinical_asl_pipeline import main_pipeline

def run_pipeline(inputdir, outputdir):
    print(f"Running pipeline for subject in: {inputdir}")
    main_pipeline.mri_diamox_umcu_clinicalasl_cvr_imager(inputdir, outputdir)
    print("Pipeline finished successfully.")

def main():
    parser = argparse.ArgumentParser(description="Run Clinical ASL CVR Pipeline")
    parser.add_argument("inputdir", type=str, help="Path to subject directory containing DICOM data")
    parser.add_argument("outputdir", type=str, help="Path to output directory for DICOM results for PACS and IMAGER platform")
    args = parser.parse_args()
    run_pipeline(args.inputdir, args.outputdir)

if __name__ == "__main__":
    main()