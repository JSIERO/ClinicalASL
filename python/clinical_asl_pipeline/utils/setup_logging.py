"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Main pipeline script for processing ASL MRI data.
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMC Utrecht (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    This code setup up the logging of the analysis
License: BSD 3-Clause License
"""
import logging
import os

def setup_logging(outputdir, log_filename='clinicalasl.log'):
    os.makedirs(outputdir, exist_ok=True)
    log_file = os.path.join(outputdir, log_filename)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, mode='w')
        ],
        force=True  # ensure it always reconfigures cleanly
    )

    logging.info(f"Logging initialized. Log file: {log_file}")