import argparse
import logging
import os
import json
from clinical_asl_pipeline import main_pipeline
from clinical_asl_pipeline.utils.load_parameters import load_parameters
from clinical_asl_pipeline.utils.setup_logging import setup_logging
from clinical_asl_pipeline.__version__ import __version__ as TOOL_VERSION
from clinical_asl_pipeline.utils.banner import log_pipeline_banner

def run_pipeline(inputdir, outputdir, inference_method=None, config_path=None):
    print(f"Running pipeline for subject in: {inputdir}")

    # Fallback default config if not supplied
    if config_path is None:
        config_path = os.path.join(os.path.dirname(__file__), 'clinical_asl_pipeline', 'config', 'config_default.json')

    # Load processing parameters (default or from user-supplied config.json)
    ANALYSIS_PARAMETERS = load_parameters(config_path=config_path)

    # Rule: If inference_method is not provided (None), use the one from config.json
    if inference_method is None:
        inference_method = ANALYSIS_PARAMETERS.get("inference_method", "basil")  # Default to 'basil' if not in config
    else:
        # Override the config's inference_method with the user-provided value
        ANALYSIS_PARAMETERS["inference_method"] = inference_method

    # Pretty print config to terminal (summary version)
    print("\n=== Key Configuration Parameters ===")
    print(f"{'Parameter':<20} | {'Value':<30}")
    print("-" * 50)
    for key, value in ANALYSIS_PARAMETERS.items():
        if isinstance(value, (str, int, float, bool)):
            print(f"{key:<20} | {str(value):<30}")
        elif isinstance(value, list) and len(value) <= 4:
            print(f"{key:<20} | {', '.join(map(str, value)):<30}")
    print("-" * 50 + "\n")
    
    # Full config to log file only
    logging.info("Full config.json contents:\n%s", 
                json.dumps(ANALYSIS_PARAMETERS, indent=4))

    # Save config copy
    with open(os.path.join(outputdir, 'config_used.json'), 'w') as f:
        json.dump(ANALYSIS_PARAMETERS, f, indent=4)

    # Run main pipeline
    main_pipeline.mri_diamox_umcu_clinicalasl_cvr(inputdir, outputdir, ANALYSIS_PARAMETERS)

    logging.info("ClinicalASL pipeline finished successfully.")
    print("Pipeline finished successfully.")

def main():
    """Main entry point for the pipeline."""
    parser = argparse.ArgumentParser(
        description="Run Clinical ASL CVR Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    Examples:
    python run_pipeline.py /input /output
    python run_pipeline.py /input /output --inference-method basil
    python run_pipeline.py /input /output --inference-method ssvb --config /path/to/config.json
        """
    )
    
    parser.add_argument("inputdir", type=str, 
                        help="Path to subject directory containing DICOM data")
    parser.add_argument("outputdir", type=str, 
                        help="Path to output directory for DICOM results for PACS and IMAGER platform")
    parser.add_argument("--inference-method", type=str, default=None, 
                        choices=["basil", "ssvb"],
                        help="Optional input to choose inference method for fitting: 'basil' or 'ssvb'. If not provided, uses the value from config.json.")
    parser.add_argument("--config", type=str, default=None,
                        help="Optional path to config.json with processing parameters")
    parser.add_argument("--version", action="version", version=f"ClinicalASL {TOOL_VERSION}")
    
    args = parser.parse_args()
    
    try:
        run_pipeline(args.inputdir, args.outputdir, 
                    inference_method=args.inference_method, 
                    config_path=args.config)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()