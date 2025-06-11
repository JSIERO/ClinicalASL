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
import subprocess
import logging
import nibabel as nib
import numpy as np
from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti
from clinical_asl_pipeline.utils.dilate_mask import dilate_mask

def run_bet_mask(m0_path, mask_output_path):
    #
    # Run HD-BET CLI on the given M0 image, save result as expected mask_path.
    #
    # Args:
    #    m0_path (str): Path to M0.nii.gz image.
    #    mask_output_path (str): Target path for brain mask, e.g. M0_brain_mask.nii.gz.
    #    output_folder (str): Folder where HD-BET will write its outputs (phase_data['ASLdir']).

    #Returns:
    #    mask (np.ndarray): Boolean brain mask.
    #    nanmask (np.ndarray): Nan-masked brain mask.
    

    # Build HD-BET CLI command
    device = 'cpu'

    cmd = f"hd-bet -i {m0_path} -o {mask_output_path} -device {device} --disable_tta --save_bet_mask"

    # Log and run command
    logging.info(f"Running brain masking with HD-BET CLI: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

    # HD-BET will create:
    # - mask_output_path_bet.nii.gz → actual mask (we want this)
    # - mask_output_path → masked image (we want to delete this)

    # Build path to HD-BET brain mask
    mask_bet_path = mask_output_path.replace('.nii.gz', '_bet.nii.gz')

    # Load brain mask
    mask = nib.load(mask_bet_path).get_fdata()
    mask = mask > 0  # Ensure binary mask



    # Save undilated mask for reference
    undil_mask_output_path = mask_output_path.replace('.nii.gz', '_undil.nii.gz')
    save_data_nifti(mask, undil_mask_output_path, m0_path, 1)
    logging.info(f"Saved undilated mask to: {undil_mask_output_path}")

    # Dilate mask
    mask_dil = dilate_mask(mask, '2D', iterations=1, conservative=False)
    nanmask = np.where(mask_dil, 1.0, np.nan)

    # Save final dilated mask
    save_data_nifti(mask_dil, mask_output_path, m0_path, 1)
    logging.info(f"Saved dilated mask to: {mask_output_path}")

    # Optionally delete unwanted HD-BET masked image
    hd_bet_masked_image = mask_output_path
    if os.path.exists(hd_bet_masked_image) and hd_bet_masked_image != mask_bet_path:
        try:
            os.remove(hd_bet_masked_image)
            logging.info(f"Deleted unwanted HD-BET masked image: {hd_bet_masked_image}")
        except Exception as e:
            logging.warning(f"Could not delete {hd_bet_masked_image}: {e}")

    return mask, nanmask