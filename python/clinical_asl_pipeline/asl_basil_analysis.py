"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL BASIL analysis module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Function for running BASIL (Bayesian Inference for Arterial Spin Labeling) analysis using FSL's oxford_asl tool.

License: BSD 3-Clause License
"""

import subprocess
import time
import logging

def asl_basil_analysis(
    subject,
    location_asl_labelcontrol_pld_nifti,
    location_m0,
    location_mask,
    output_map,
    pld_list,
    artoff=None,
    spatialoff=None
):
# This function runs the BASIL (Bayesian Inference for Arterial Spin Labeling) analysis
# using the FSL oxford_asl tool. It prepares the command line arguments based on the
# provided subject parameters and input file locations, then executes the analysis.
# Parameters:
# subject: dict containing subject information including:
#   - 'T1t': T1 tissue relaxation time in seconds
#   - 'T1b': T1 blood relaxation time in seconds
#   - 'tau': bolus duration in seconds
#   - 'TR_M0': list or tuple, first element is TR for M0 scan in seconds
#   - 'alpha': labeling efficiency (unitless)
#   - 'slicetime': slice timing in milliseconds
# location_asl_labelcontrol_pld_nifti: str, path to ASL label/control PLD NIfTI file
# location_m0: str, path to M0 image file
# location_mask: str, path to brain mask file
# output_map: str, path for output directory
# pld_list: list of PLDs (post-labeling delays) in seconds
# artoff: optional, set to "artoff" to disable arterial component modeling
# spatialoff: optional, set to "spatialoff" to disable spatial regularization

    # Generate comma-separated PLD string
    pld_string = ",".join([f"{pld:.5g}" for pld in pld_list])

    # Arterial component off (optional)
    artoff_string = " --artoff" if artoff == "artoff" else ""

    # Spatial regularization (optional)
    spatial_string = "off" if spatialoff == "spatialoff" else "on"

    # Extract parameter values and convert to string
    T1t = str(subject['T1t'])
    T1b = str(subject['T1b'])
    tau = str(subject['tau'])
    TR_M0 = str(subject['TR_M0'][0])
    alpha = str(subject['alpha'])
    slicetime = str(subject['slicetime'] / 1000)  # convert ms to seconds

    # Timing the execution
    start_time = time.time()

    # Build oxford_asl command
    cmd = (
        f"oxford_asl "
        f"-i {location_asl_labelcontrol_pld_nifti} "
        f"-c {location_m0} "
        f"-m {location_mask} "
        f"-o {output_map} "
        f"{artoff_string} "
        f"--spatial={spatial_string} "
        f"--bolus={tau} "
        f"--slicedt={slicetime} "
        f"--t1={T1t} "
        f"--t1b={T1b} "
        f"--t1t={T1t} "
        f"--plds={pld_string} "
        f"--tr={TR_M0} "
        f"--alpha={alpha} "
        f"--iaf=ct "
        f"--ibf=tis "
        f"--casl "
        f"--fixbolus "
        f"--cmethod voxel "
        f"--cgain 1.00"
    )

    # Run command
    logging.info("Running BASIL analysis...")
    subprocess.run(cmd, shell=True, check=True)

    logging.info("BASIL analysis finished")
    elapsed = round(time.time() - start_time, 2)
    logging.info(f"..this took: {elapsed} s")
