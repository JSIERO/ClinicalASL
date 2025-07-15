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
import sys
import os
import json
from clinical_asl_pipeline import main_pipeline
from clinical_asl_pipeline.utils.load_parameters import load_parameters
from clinical_asl_pipeline.utils.setup_logging import setup_logging
from clinical_asl_pipeline.__version__ import __version__ as TOOL_VERSION
from clinical_asl_pipeline.utils.banner import log_pipeline_banner

def format_config_for_display(config_dict):
    """Format configuration parameters in a readable table format."""
    lines = []
    lines.append("=== Key Configuration Parameters ===")
    lines.append(f"{'Parameter':<20} | {'Value':<30}")
    lines.append("-" * 50)
    
    for key, value in config_dict.items():
        if isinstance(value, (str, int, float, bool)):
            lines.append(f"{key:<20} | {str(value):<30}")
        elif isinstance(value, list) and len(value) <= 4:
            lines.append(f"{key:<20} | {', '.join(map(str, value)):<30}")
    
    lines.append("-" * 50)
    return "\n".join(lines)

def run_pipeline(inputdir, outputdir, workingdir, inference_method=None, config_path=None):
    print(f"Running pipeline for subject in: {inputdir}")

    # Fallback default config if not supplied
    if config_path is None:
        config_path = os.path.join(os.path.dirname(__file__), 'clinical_asl_pipeline', 'config', 'config_default.json')

    # Load processing parameters (default or from user-supplied config.json)
    ANALYSIS_PARAMETERS = load_parameters(config_path=config_path)

    # Rule: If inference_method is not provided (None), use the one from config.json
    if inference_method is None:
        inference_method = ANALYSIS_PARAMETERS.get("inference_method", "ssvb")  # Default to 'ssvb' if not in config
    else:
        # Override the config's inference_method with the user-provided value
        ANALYSIS_PARAMETERS["inference_method"] = inference_method

    # Setup logging (also save to output dir)
    setup_logging(outputdir)  
    log_pipeline_banner(TOOL_VERSION)
    logging.info(f"ClinicalASL pipeline started.")
    logging.info(f"Input: {inputdir}")
    logging.info(f"Output: {outputdir}")
    logging.info(f"Working Directory: {workingdir}")
    logging.info(f"Config used: {config_path}")
    logging.info(f"Inference method: {inference_method}")  # Log the chosen method

    # Format config for display and logging
    formatted_config = format_config_for_display(ANALYSIS_PARAMETERS)
    logging.info(f"Configuration parameters:\n{formatted_config}")
    
    # Save config.json copy
    with open(os.path.join(outputdir, 'config_used.json'), 'w') as f:
        json.dump(ANALYSIS_PARAMETERS, f, indent=4)

    # Run main pipeline
    main_pipeline.mri_diamox_umcu_clinicalasl_cvr(inputdir, outputdir, workingdir, ANALYSIS_PARAMETERS)

    logging.info("ClinicalASL pipeline finished successfully.")
    print("Pipeline finished successfully.")

def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description="Run Clinical ASL CVR Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    Examples:
    python run_pipeline.py /input /output /working
    python run_pipeline.py /input /output /working --inference-method ssvb --config /path/to/config.json
    python run_pipeline.py /input /output /working --inference-method vaby for BASIL-like output
    """
    )
    
    parser.add_argument("inputdir", type=str, 
                        help="Path to subject directory containing DICOM data")
    parser.add_argument("outputdir", type=str, 
                        help="Path to output directory for DICOM results for PACS and IMAGR platform") 
    parser.add_argument("workingdir", type=str, 
                        help="Path to working directory for intermediate files (e.g., /tmp/workingdir")
    parser.add_argument("--inference-method", type=str, 
                        default='ssvb', 
                        choices=["ssvb", "vaby"],
                        help="Optional input to choose inference method for fitting: 'ssvb' or 'vaby'. If not provided, uses the value from config.json.")
    parser.add_argument("--config", type=str, default=None,
                        help="Optional path to config.json with processing parameters")
    parser.add_argument("--version", action="version", version=f"ClinicalASL {TOOL_VERSION}")
    
    args = parser.parse_args()
    
    try:
        run_pipeline(args.inputdir, args.outputdir, args.workingdir,
                    inference_method=args.inference_method, 
                    config_path=args.config)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()