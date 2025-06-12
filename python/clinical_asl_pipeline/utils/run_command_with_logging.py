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

    # Start subprocess
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # Read stdout and stderr line by line
    while True:
        # Read one line from stdout
        stdout_line = process.stdout.readline()
        if stdout_line:
            stdout_line = stdout_line.rstrip()
            print(stdout_line)  # print live to terminal
            logging.info(stdout_line)

        # Read one line from stderr
        stderr_line = process.stderr.readline()
        if stderr_line:
            stderr_line = stderr_line.rstrip()
            print(stderr_line)  # print live to terminal
            logging.warning(stderr_line)

        # Check if process is finished
        if process.poll() is not None:
            break

    # Check return code
    retcode = process.wait()
    if retcode != 0:
        raise subprocess.CalledProcessError(retcode, cmd)

    logging.info(f"Command finished with return code {retcode}")

