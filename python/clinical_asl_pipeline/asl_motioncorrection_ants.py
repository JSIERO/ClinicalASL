"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL data preparation module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Perform motion correction on 4D datasets using a reference image.
    Method is by ANTs

License: BSD 3-Clause License
"""

import os
import ants 
import logging

def asl_motioncorrection_ants(inputdata, refdata, outputdata):
    # Perform motion correction using ANTs
    # inputdata: path to the input ASL data in NIfTI format
    # refdata: path to the reference image for motion correction
    # outputdata: path to save the motion-corrected output data in NIfTI format
    logging.info(f"Perform motion correction (using ANTs): input: {os.path.basenam(inputdata)}  reference: {os.path.basenam(refdata)} ")
    results_dict = ants.motion_correction(
        ants.image_read(inputdata),
        ants.image_read(refdata),
        type_of_transform="DenseRigid",
        aff_metric="mattes",
        smoothing_in_mm=True,
        singleprecision=True,
        verbose=True
    )
    ants.image_write(results_dict["motion_corrected"], outputdata)
