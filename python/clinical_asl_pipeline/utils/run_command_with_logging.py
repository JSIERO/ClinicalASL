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
import threading

import subprocess
import logging
import threading

def run_command_with_logging(cmd):
    logging.info(f"Running command: {cmd}")

    # Define known harmless phrases to suppress (can expand this list)
    suppress_phrases = [
        'dcmodify: Modify DICOM files',
        'usage: dcmodify',
    ]

    def stream_output(stream, log_function, log_in_file=True):
        for line in iter(stream.readline, ''):
            if not line:
                break
            line = line.rstrip()
            print(line)

            # Suppress harmless messages
            if any(phrase in line for phrase in suppress_phrases):
                continue

            if log_in_file:
                log_function(line)

    # stdout direct to terminal → progress bar works
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=None,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    # Start thread to read stderr
    stderr_thread = threading.Thread(target=stream_output, args=(process.stderr, logging.info, False))
    stderr_thread.start()

    # Wait for command to finish
    retcode = process.wait()
    stderr_thread.join()

    if retcode != 0:
        raise subprocess.CalledProcessError(retcode, cmd)

    logging.info(f"Command finished with return code {retcode}")
