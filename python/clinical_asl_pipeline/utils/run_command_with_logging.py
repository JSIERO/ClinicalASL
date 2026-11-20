"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Run an external command (dcm2niix, hd-bet, qasl) with logging.

Repository: https://github.com/JSIERO/ClinicalASL
Author: Jeroen Siero, UMCU (University Medical Center Utrecht), The Netherlands
License: BSD 3-Clause License
"""
import subprocess
import logging
import threading
from collections import deque

STDERR_TAIL_LINES = 40  # stderr lines written to the log file when a command fails


def run_command_with_logging(cmd):
    # stdout goes straight to the console (progress bars work); stderr is echoed to the console and its last
    # STDERR_TAIL_LINES lines are kept, to be logged (console + clinicalasl.log) only if the command fails.
    logging.info(f"Running command: {cmd}")
    tail = deque(maxlen=STDERR_TAIL_LINES)

    def stream_stderr(stream):
        for line in iter(stream.readline, ''):
            line = line.rstrip()
            print(line)
            tail.append(line)

    process = subprocess.Popen(cmd, shell=True, stdout=None, stderr=subprocess.PIPE, text=True, bufsize=1)
    stderr_thread = threading.Thread(target=stream_stderr, args=(process.stderr,))
    stderr_thread.start()
    retcode = process.wait()
    stderr_thread.join()

    if retcode != 0:
        logging.error(f"Command failed with return code {retcode}; last {len(tail)} stderr line(s):\n" + "\n".join(tail))
        raise subprocess.CalledProcessError(retcode, cmd)
    logging.info(f"Command finished with return code {retcode}")
