"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Main pipeline script for processing ASL MRI data.
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMC Utrecht (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    This code loads in the study or protocol specific scan and analysis parameters
License: BSD 3-Clause License
"""
import logging
import os
import json

def load_parameters(config_path=None):
    # If no config_path provided, use default
    if not config_path:
        config_path = os.path.join(os.path.dirname(__file__), 'config', 'config_default.json')

    # Now safe to check if exists
    if not os.path.exists(config_path):
        logging.warning(f"Config file '{config_path}' not found. Using hardcoded defaults.")
        hardcoded_defaults =  {
        "version": "v1.0-MRI_DIAMOX_MDLL_preACZ_postACZ-2025",
        "ASL scan": "multi-delay Look-Locker",
        "tau": 2,
        "N_BS": 4,
        "readout": "2D",
        "labeleff": 0.85,
        "lambda": 0.9,
        "T1t": 1.3,
        "T1b": 1.65,
        "FWHM": 6,
        "outlier_factor": 2.5,
        "range_cbf": [0, 100],
        "range_cvr": [-50, 50],
        "range_AAT": [0, 3],
        "range_ATA": [0, 125],
        "inference_method": "ssvb",
        "device": "cpu",
        "ASL_CONTEXT": ["baseline", "stimulus"],
        "context_study_tags": ["preACZ", "postACZ"],
        "dicomseries_description_patterns": ["*SOURCE*ASL*", "SWIP*ASL*"]
        }
        
        logging.info(f"Config parameters (hardcoded):\n{json.dumps(hardcoded_defaults, indent=2)}")
        return hardcoded_defaults
    else:
        with open(config_path, 'r') as f:
            params = json.load(f)
        version = params.get('version', 'unknown')
        logging.info(f"Loaded processing parameters from '{config_path}' (version: {version})")
        # Log the full config content:
        logging.info(f"Config parameters:\n{json.dumps(params, indent=2)}")
        return params