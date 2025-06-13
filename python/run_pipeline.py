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
import logging
import os
from clinical_asl_pipeline import main_pipeline
from clinical_asl_pipeline.utils.load_parameters import load_parameters
from clinical_asl_pipeline.utils.setup_logging import setup_logging
from clinical_asl_pipeline.__version__ import __version__ as TOOL_VERSION
from clinical_asl_pipeline.utils.banner import log_pipeline_banner

def run_pipeline(inputdir, outputdir, config_path=None):
    print(f"Running pipeline for subject in: {inputdir}")

    # Fallback default config if not supplied
    if config_path is None:
        config_path = os.path.join(os.path.dirname(__file__), 'clinical_asl_pipeline', 'config', 'config_default.json')

    # Load processing parameters (default or from user-supplied config.json)
    ANALYSIS_PARAMETERS = load_parameters(config_path=config_path)

    # Setup logging (also save to output dir)
    setup_logging(outputdir)  
    log_pipeline_banner(TOOL_VERSION)
    logging.info(f"Pipeline config version: {ANALYSIS_PARAMETERS.get('version', 'unknown')}")
    logging.info(f"ClinicalASL pipeline started.")
    logging.info(f"Input: {inputdir}")
    logging.info(f"Output: {outputdir}")
    logging.info(f"Config used: {config_path}")
    logging.info(f"Pipeline config version: {ANALYSIS_PARAMETERS.get('version', 'unknown')}")

    # Run main pipeline
    main_pipeline.mri_diamox_umcu_clinicalasl_cvr(inputdir, outputdir, ANALYSIS_PARAMETERS)

    logging.info("ClinicalASL pipeline finished successfully.")
    print("Pipeline finished successfully.")

def main():
    parser = argparse.ArgumentParser(description="Run Clinical ASL CVR Pipeline")
    parser.add_argument("inputdir", type=str, help="Path to subject directory containing DICOM data")
    parser.add_argument("outputdir", type=str, help="Path to output directory for DICOM results for PACS and IMAGER platform")
    parser.add_argument("--config", type=str, default=None, help="Optional path to config.json with processing parameters")
    args = parser.parse_args()

    run_pipeline(args.inputdir, args.outputdir, config_path=args.config)

if __name__ == "__main__":
    main()