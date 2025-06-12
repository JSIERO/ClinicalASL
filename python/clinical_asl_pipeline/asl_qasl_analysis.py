"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline
Repository: https://github.com/JSIERO/ClinicalASL

Script Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Algorithm code author: Thomas Kirk, Quantified Imaging

Description:
    Performs QASL - by Quantified Imaging (T.Kirk) (Quantitative Arterial Spin Labeling) analysis on ASL data using the Oxford ASL toolbox.
    This script provides the asl_qasl_analysis() function, which builds and runs a command-line call to the QASL tool,
    passing all relevant parameters for quantification.

License: BSD 3-Clause License
"""

import time
import logging
from clinical_asl_pipeline.utils.run_command_with_logging import run_command_with_logging

def asl_qasl_analysis(
    subject,
    ANALYSIS_PARAMETERS, 
    location_asl_labelcontrol_pld_nifti,
    location_m0,
    location_mask,
    output_map,
    pld_list,
    inference_method='ssvb', # or basil, for BASIL method output
    artoff=None,
):
    # Perform QASL analysis on ASL data using the Oxford ASL toolbox.
    # Parameters:
    # subject: dict containing subject information including, can be different per context tag   
    #   - 'T1t': T1 tissue relaxation time in seconds
    #   - 'T1b': T1 blood relaxation time in seconds
    #   - 'tau': bolus duration in seconds
    #   - 'TR_M0': list or tuple, first element is TR for M0 scan
    #   - 'alpha': labeling efficiency
    #   - 'slicetime': slice timing in milliseconds
    # location_asl_labelcontrol_pld_nifti: path to ASL label/control NIfTI file
    # location_m0: path to M0 NIfTI file
    # location_mask: path to brain mask NIfTI file
    # output_map: output directory for QASL results
    # pld_list: list of post-labeling delays (PLDs) in seconds
    # inference_method: inference method for QASL ('ssvb' or 'basil')
    # artoff: optional, set to 'artoff' to disable arterial component modeling
    #
    # This function builds and runs a command-line call to the QASL tool,
    # passing all relevant parameters for quantification. It times the execution,
    # prints progress messages, and ensures the command is run with error checking.    

    # Generate comma-separated PLD string
    pld_string = ",".join([f"{pld:.5g}" for pld in pld_list])

    # Arterial component off (optional)
    artoff_string = " --artoff" if artoff == "artoff" else ""

    # Extract parameter values and convert to string
    T1t = str(ANALYSIS_PARAMETERS['T1t'])
    T1b = str(ANALYSIS_PARAMETERS['T1b'])
    tau = str(ANALYSIS_PARAMETERS['tau'])
    readout = str(ANALYSIS_PARAMETERS['readout']) # 2D or 3D
    alpha = str(subject['alpha'])
    TR_m0 = str(subject['TR_M0'][0])
    slicetime = str(subject['slicetime'] / 1000)  # convert ms to seconds
    # Timing the execution
    start_time = time.time()

    # Build qasl command
    cmd = (
        f"qasl -i {location_asl_labelcontrol_pld_nifti} "
        f"-c {location_m0} "
        f"-m {location_mask} "
        f"-o {output_map} "
        f" --inference-method={inference_method} "
        f"{artoff_string} "
        f"--bolus={tau} "
        f"--slicedt={slicetime} "
        f"--t1={T1t} "
        f"--t1b={T1b} "
        f"--t1t={T1t} "
        f"--plds={pld_string} "
        f"--tr={TR_m0} "
        f" --alpha={alpha} "
        f"--iaf=ct "
        f"--ibf=tis "
        f"--casl "
        f"--cgain 1.00 "
        f"--readout={readout} "
        f"--save-calib "
        f"--overwrite"
    )

    # Run command
    logging.info("Running QASL analysis...")
    run_command_with_logging(cmd)

    logging.info("QASL analysis finished")
    elapsed = round(time.time() - start_time, 2)
    logging.info(f"..this took: {elapsed} s")
