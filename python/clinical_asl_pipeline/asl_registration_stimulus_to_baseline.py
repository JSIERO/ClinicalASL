"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

DICOM to NIfTI conversion module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Function for registering ASL stimulus data to baseline data using ANTsPy.

License: BSD 3-Clause License
"""

import os
import logging
import shutil
import ants

def asl_registration_stimulus_to_baseline(subject):
    # Register post-ACZ ASL data to pre-ACZ ASL data using ANTsPy
    # subject is a dictionary containing paths to the necessary files.
    # The dictionary should contain the following keys:
    # 'baseline': {
    #     'M0_path', 'CBF_path', 'AAT_path', 'ATA_path', 'mask_path'
    # },
    # 'stimulus': {
    #     'M0_path', 'CBF_path', 'AAT_path', 'ATA_path', 'mask_path',
    #     'M0_2baseline_path', 'CBF_2baseline_path', 'AAT_2baseline_path',
    #     'ATA_2baseline_path', 'mask_2baseline_path'
    # },
    # 'ASLdir'
    # resulting transform will be saved in 'ASLdir' as 'rigid_stimulus_to_baseline.mat'

    # Load fixed and moving images for registration
    logging.info("Registration M0 stimulus to baseline data (ANTsPy) *********************************************************************")

    fixed = ants.image_read(subject['baseline']['M0_path'])
    moving = ants.image_read(subject['stimulus']['M0_path'])
    
    # Run registration
    reg = ants.registration(fixed=fixed, moving=moving, type_of_transform='Rigid', metric='Mattes', reg_iterations=(1000, 500, 250, 100)) 
    interpolator = 'bSpline'
    # Save transformed moving image
    ants.image_write(reg['warpedmovout'], subject['stimulus']['M0_2baseline_path'])

    # Apply same transform to CBF, AAT, ATA, mask
    def apply_transform(moving_path, reference_path, output_path, transformlist, interpolation):
        moving_img = ants.image_read(moving_path)
        reference_img = ants.image_read(reference_path)
        warped = ants.apply_transforms(fixed=reference_img, moving=moving_img,
                                        transformlist=transformlist,
                                        interpolator=interpolation)
        ants.image_write(warped, output_path)

    # CBF
    logging.info("Registration CBF stimulus to baseline (ANTsPy)")
    apply_transform(subject['stimulus']['QASL_CBF_path'], subject['baseline']['QASL_CBF_path'], subject['stimulus']['CBF_2baseline_path'], reg['fwdtransforms'], interpolator)

    # AAT
    logging.info("Registration AAT stimulus to baseline (ANTsPy)")
    apply_transform(subject['stimulus']['QASL_AAT_path'], subject['baseline']['QASL_AAT_path'], subject['stimulus']['AAT_2baseline_path'], reg['fwdtransforms'], interpolator)

    # ATA
    if subject['dicom_typetags_by_context']['baseline'].__contains__('ATA'): # only if ATA is present in baseline, e.g. ATA not yet generated for vTR data
        logging.info("Registration ATA stimulus to baseline (ANTsPy)")
        apply_transform(subject['stimulus']['QASL_ATA_path'], subject['baseline']['QASL_ATA_path'], subject['stimulus']['ATA_2baseline_path'], reg['fwdtransforms'], interpolator)

    # Mask (NearestNeighbor interpolation)
    logging.info("Registration mask stimulus to baseline (ANTsPy)")
    apply_transform(subject['stimulus']['mask_path'], subject['baseline']['mask_path'], subject['stimulus']['mask_2baseline_path'], reg['fwdtransforms'], 'nearestNeighbor')

    transform_path = reg['fwdtransforms'][0]

    # Choose your save location for transform (ITK format)
    output_transform_path = os.path.join(subject['ASLdir'], 'rigid_stimulus_to_baseline.mat')

    # Copy it to save location
    shutil.copy(transform_path, output_transform_path)