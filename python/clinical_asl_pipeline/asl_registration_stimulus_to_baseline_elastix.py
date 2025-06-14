"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Performs registration post-ACZ ASL data to pre-ACZ ASL data using Elastix

License: BSD 3-Clause License
"""

import os
import importlib.resources as pkg_resources
import tempfile
import shutil
import warnings
import logging
from clinical_asl_pipeline.utils.run_command_with_logging import run_command_with_logging

def asl_registration_prepostACZ_elastix(subject):
    # Function to register post-ACZ ASL data to pre-ACZ ASL data using Elastix
    # subject is a dictionary containing paths to the necessary files.
    # The dictionary should contain the following
    # keys: 
    # 'preACZ_T1fromM0_path', 'postACZ_T1fromM0_path', 'preACZ_mask_path',
    # 'postACZ_CBF_path', 'preACZ_CBF_path', 'postACZ_AAT_path',
    # 'preACZ_AAT_path', 'postACZ_ATA_path', 'preACZ_ATA_path',
    # 'postACZ_T1fromM0_2preACZ_path', 'postACZ_CBF_2preACZ_path',
    # 'postACZ_AAT_2preACZ_path', 'postACZ_ATA_2preACZ_path',
    # 'postACZ_mask_2preACZ_path', 'ASLdir', 'elastix_parameter_file'
    
    try:
        # Get the resource path (for 3.9+)
        elastix_file_path = pkg_resources.files('clinical_asl_pipeline.elastixfiles').joinpath('Par0001rigid_6DOF_MI_NIFTIGZ.txt')  # Elastix registration parameter file

        # Copy to temp location to get a real filesystem path
        with tempfile.NamedTemporaryFile(delete=False, suffix='.txt') as tmpfile:
            shutil.copy(elastix_file_path, tmpfile.name)
            subject['elastix_parameter_file'] = tmpfile.name

    except FileNotFoundError:
        warnings.warn(f"Missing Elastix file: {subject['elastix_parameter_file']}")
        return
    
    # Elastix registration
    logging.info("Registration T1fromM0 postACZ to preACZ data *********************************************************************")
    run_command_with_logging(f"elastix -f {subject['preACZ']['T1fromM0_path']} -m {subject['postACZ']['T1fromM0_path']} -fMask {subject['preACZ']['mask_path']} -p {subject['elastix_parameter_file']} -loglevel error -out {subject['ASLdir']}")
    os.replace(os.path.join(subject['ASLdir'], 'result.0.nii.gz'), subject['postACZ']['T1fromM0_2preACZ_path'])
    run_command_with_logging(f"fslcpgeom {subject['preACZ']['T1fromM0_path']} {subject['postACZ']['T1fromM0_2preACZ_path']} -d")

    logging.info("Registration CBF postACZ to preACZ *********************************************************************")
    run_command_with_logging(f"transformix -in {subject['postACZ']['QASL_CBF_path']} -out {subject['ASLdir']} -tp {os.path.join(subject['ASLdir'], 'TransformParameters.0.txt')} -loglevel error")
    os.replace(os.path.join(subject['ASLdir'], 'result.nii.gz'), subject['postACZ']['CBF_2preACZ_path'])
    run_command_with_logging(f"fslcpgeom {subject['preACZ']['QASL_CBF_path']} {subject['postACZ']['CBF_2preACZ_path']} -d")

    logging.info("Registration AAT postACZ to preACZ *********************************************************************")
    run_command_with_logging(f"transformix -in {subject['postACZ']['QASL_AAT_path']} -out {subject['ASLdir']} -tp {os.path.join(subject['ASLdir'], 'TransformParameters.0.txt')} -loglevel error")
    os.replace(os.path.join(subject['ASLdir'], 'result.nii.gz'), subject['postACZ']['AAT_2preACZ_path'])
    run_command_with_logging(f"fslcpgeom {subject['preACZ']['QASL_AAT_path']} {subject['postACZ']['AAT_2preACZ_path']} -d")
    
    logging.info("Registration ATA postACZ to preACZ *********************************************************************")
    run_command_with_logging(f"transformix -in {subject['postACZ']['QASL_ATA_path']} -out {subject['ASLdir']} -tp {os.path.join(subject['ASLdir'], 'TransformParameters.0.txt')} -loglevel error")
    os.replace(os.path.join(subject['ASLdir'], 'result.nii.gz'), subject['postACZ']['ATA_2preACZ_path'])
    run_command_with_logging(f"fslcpgeom {subject['preACZ']['QASL_ATA_path']} {subject['postACZ']['ATA_2preACZ_path']} -d")

    logging.info("Registration mask postACZ to preACZ *********************************************************************")
    with open(os.path.join(subject['ASLdir'], 'TransformParameters.0.txt'), 'r') as f:
        lines = f.readlines()
    with open(os.path.join(subject['ASLdir'], 'TransformParameters.0.NN.txt'), 'w') as f:
        for line in lines:
            if 'FinalBSplineInterpolator' in line:
                f.write('(ResampleInterpolator "FinalNearestNeighborInterpolator")\n')
            else:
                f.write(line)

    run_command_with_logging(f"transformix -in {subject['postACZ']['mask_path']} -out {subject['ASLdir']} -tp {os.path.join(subject['ASLdir'], 'TransformParameters.0.NN.txt')} -loglevel error")
    os.replace(os.path.join(subject['ASLdir'], 'result.nii.gz'), subject['postACZ']['mask_2preACZ_path'])
    run_command_with_logging(f"fslcpgeom {subject['preACZ']['mask_path']} {subject['postACZ']['mask_2preACZ_path']} -d")
    
    os.replace(subject['elastix_parameter_file']) # clean up temporary elastix parameter file

