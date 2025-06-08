#!/usr/bin/env python3
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