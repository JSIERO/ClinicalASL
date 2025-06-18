"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL T1-from-M0 computation module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Functions for brain extraction using HD-BET CLI.

License: BSD 3-Clause License
"""
import os
import logging
import nibabel as nib
import numpy as np
from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti
from clinical_asl_pipeline.utils.dilate_mask import dilate_mask
from clinical_asl_pipeline.utils.run_command_with_logging import run_command_with_logging

def run_bet_mask(inputdata_path, mask_output_path, device='cpu', extradata_path=None):
    #
    # Run HD-BET CLI on the given M0 image, save result as expected mask_path.
    #
    # Args:
    #    inputdata_path (str): Path to M0.nii.gz image.
    #    mask_output_path (str): Target path for brain mask, e.g. M0_brain_mask.nii.gz.

    #Returns:
    #    mask (np.ndarray): Boolean brain mask.
    #    nanmask (np.ndarray): Nan-masked brain mask.

    templateNII_path = inputdata_path
    logging.info(f"Running brain masking with HD-BET CLI:")

    # Build HD-BET CLI command
    cmd = f"hd-bet -i {inputdata_path} -o {mask_output_path} -device {device} --disable_tta --save_bet_mask"

    if extradata_path and os.path.exists(extradata_path):
        # Load input and extra data
        inputdata = nib.load(inputdata_path).get_fdata()
        extradata = nib.load(extradata_path).get_fdata()

        # Ensure inputdata is 4D for concatenation
        if inputdata.ndim == 3:
            inputdata = inputdata[..., np.newaxis]

        # Concatenate and sum across time dimension
        combineddata = np.concatenate((inputdata, extradata), axis=3).sum(axis=3)

        # Save temporary image for brain extraction
        temp_bet_path = os.path.join(os.path.dirname(mask_output_path), "temp_for_bet.nii.gz")
        save_data_nifti(combineddata, temp_bet_path, templateNII_path, 1)
        logging.info(f"Using combined data set for brain masking: sum of  {inputdata_path} and {extradata_path}")
        inputdata_path = temp_bet_path
        # run command HD-BET
        run_command_with_logging(cmd)
        # remove temp file of combined data for full covering brain mask
        os.remove(temp_bet_path) 
    else:
        run_command_with_logging(cmd)

    # HD-BET will create:
    # - mask_output_path_bet.nii.gz → actual mask (we want this)
    # - mask_output_path → masked image (we want to delete this)

    # Build path to HD-BET brain mask
    mask_bet_path = mask_output_path.replace('.nii.gz', '_bet.nii.gz')

    # Load brain mask
    mask = nib.load(mask_bet_path).get_fdata()
    mask = mask > 0  # Ensure binary mask

    # Dilate mask 1 voxels: using 3D and conservative mode for a tight mask
    mask_dil = dilate_mask(mask, '3D', iterations=1, conservative=True)

    nanmask = np.where(mask_dil, 1.0, np.nan)

    # Save final dilated mask
    save_data_nifti(mask_dil, mask_output_path, templateNII_path, 1)
    logging.info(f"Saved dilated mask to: {mask_output_path}")

    # Delete unwanted HD-BET masked image
    os.remove(mask_bet_path)
    logging.info(f"Deleted unwanted HD-BET masked image: {mask_bet_path}")

    return mask, nanmask