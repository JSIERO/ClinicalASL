"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Main pipeline script for processing ASL MRI data.
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMC Utrecht (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    This utiliy also logs command line output by subprocess

License: BSD 3-Clause License
"""

import subprocess
import logging

def run_command_with_logging(cmd):
    logging.info(f"Running command: {cmd}")

    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # Log stdout lines if any
    if result.stdout:
        for line in result.stdout.strip().splitlines():
            logging.info(f"    {line}")

    # Log stderr lines if any
    if result.stderr:
        for line in result.stderr.strip().splitlines():
            logging.warning(f"    {line}")

    return result
