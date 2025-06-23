"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Main pipeline script for processing ASL MRI data.
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Outlier removal by rejecting volumes from the data,
    An outlier volume is one whose mean deltaM (in mask) deviates more than mean deltaM over volumes + outlier_factor x temporal std brain deltaM 
    precisely the median across the volumes is used, and mad (median absolute deviation)

License: BSD 3-Clause License
"""

import numpy as np
import logging
import nibabel as nib
from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti
from clinical_asl_pipeline.utils.append_filename import append_or


def mad(data, axis=None):
    # Median Absolute Deviation: a robust version of standard deviation."""
    return np.median(np.abs(data - np.median(data, axis=axis)), axis=axis)




def asl_outlier_removal(subject, context_tag, usermask=None):
    # Remove outlier volumes from ASL deltaM signal based on mean deltaM in a combined mask.   
    # method: outlier volume is one whose mean deltaM (in mask) deviates more than mean deltaM over volumes + outlier_factor x temporal std brain deltaM 
    #         precisely the median across the volumes is used, and mad (median absolute deviation)
    # Parameters:
    #   subject: dict containing subject information including paths and parameters
    #   context_tag: string, e.g. 'baseline' or 'stimulus' context_tag for the keys in the subject dictionary to store results
    #   usermask = user supplied mask to be combined with the brainmask (standard mask)
    # Returns:
    #   outlier-removed PLD ordered NIFTIs, filename appended with '_or' before '.nii.gz'

    # Use a shorter alias for subject[context_tag]
    context_data = subject[context_tag]
    outlier_factor = subject['outlier_factor']    

    M0ASL_allPLD6D = context_data['M0ASL_allPLD']
    brainmask_path = context_data['mask_path']
    ASL_allPLD4D_path = context_data['PLDall_controllabel_path']
    nifti_template_path = context_data['templateNIFTI_path']
    NREPEATS = context_data['NREPEATS']
    NPLDS = context_data['NPLDS']
    TR = context_data['TR'] 

    # Load mask and combine with user-supplied mask if provided
    brainmask = nib.load(brainmask_path).get_fdata() > 0
    mask = brainmask if usermask is None else np.logical_and(brainmask, usermask)

    # Fetch ASL data and exclude first volume (M0)
    # shape is (x, y, z, NREPEATS, NPLDS, control/label)
    ASL_allPLD = M0ASL_allPLD6D[:, :, :, 1:, :, :]

    # Compute deltaM = label - control across time
    deltaM5D = ASL_allPLD[..., 1] - ASL_allPLD[..., 0]  # shape: (x, y, z, NREPEATS, NPLDS)
    deltaM4D = deltaM5D.sum(axis=4)  # shape: (x, y, z, NREPEATS)

    # Masked view of data: shape (N_voxels, NREPEATS)
    voxels_in_mask = deltaM4D[mask]
    voxels_in_mask = np.where(np.isinf(voxels_in_mask), 0, voxels_in_mask)

    # Mean deltaM per volume
    mean_per_volume = np.nanmean(voxels_in_mask, axis=0)

    if np.isnan(mean_per_volume).any():
        raise ValueError(f"NaNs found in deltaM means at volumes: {np.where(np.isnan(mean_per_volume))[0]}")

    ## Find outlier volumes using MAD (robust)
    mu = np.median(mean_per_volume)
    std = 1.4826 * mad(mean_per_volume)
    keep_volumes = np.abs(mean_per_volume - mu) < outlier_factor * std

    outlier_indices = np.where(~keep_volumes)[0]
    keep_indices = np.where(keep_volumes)[0]
    NREPEATS_kept = len(keep_indices)
    n_vols_per_pld = 2 * NREPEATS_kept

    context_data['outlier_volumes_indices'] = outlier_indices

    logging.info(f"===================================================================")

    # Remove outlier volumes from ASL data (4D) needed for QASL
    if outlier_indices.size > 0:
        logging.info(f"Outlier removal: Volumes removed (1-based): {(outlier_indices + 1).tolist()}")
        # Load the 4D NIfTI
        PLDall = nib.load(ASL_allPLD4D_path).get_fdata()
        x, y, z, t = PLDall.shape
        assert t == 2 * NREPEATS * NPLDS, "Time dimension doesn't match expected size."

        # Step 1: Reshape to (x, y, z, NPLDS, NREPEATS, control/label)
        PLDall_reshaped = PLDall.reshape(x, y, z, NPLDS, NREPEATS, 2)

        # Step 2: Remove unwanted repeats
        # keep_indices: list or array of repeat indices to keep (length = NREPEATS_kept)
        PLDall_kept = PLDall_reshaped[:, :, :, :, keep_indices, :]  # (x, y, z, NPLDS, NREPEATS_kept, control/label)

        # Step 3: Reshape back to 4D for saving
        # First: transpose to match original time order: PLD (outer), repeat, control/label (inner)
        # From shape: (x, y, z, NPLDS, NREPEATS_kept, 2) → (x, y, z, control/label * NREPEATS_kept * NPLDS)
        PLDall_or = np.transpose(PLDall_kept, (0, 1, 2, 3, 4, 5)).reshape(x, y, z, n_vols_per_pld * NPLDS)

        # Derive subsets
        PLD1to2_or = PLDall_or[:, :, :, :2 * n_vols_per_pld]      # first 2 PLDs, for ATA
        PLD2tolast_or = PLDall_or[:, :, :, n_vols_per_pld:]       # from PLD2 onward, for CBF

        # Save all outputs 
        logging.info(f"Saving ASL outlier-removed data: 'all PLDs for AAT'")
        logging.info(f"Saving ASL outlier-removed data: '2-to-last PLDs for CBF'")
        logging.info(f"Saving ASL outlier-removed data: '1-to-2 PLDs for ATA'")

        # update path to outlier removed corrected data, appeding '_or' to filename using append_or
        context_data['PLDall_controllabel_path'] =  append_or(context_data['PLDall_controllabel_path'])
        context_data['PLD2tolast_controllabel_path'] =  append_or(context_data['PLD2tolast_controllabel_path']) 
        context_data['PLD1to2_controllabel_path'] =  append_or(context_data['PLD1to2_controllabel_path'])

        save_data_nifti(PLDall_or, context_data['PLDall_controllabel_path'], nifti_template_path,  1, None, TR)
        save_data_nifti(PLD2tolast_or, context_data['PLD2tolast_controllabel_path'], nifti_template_path,  1, None, TR)        
        save_data_nifti(PLD1to2_or, context_data['PLD1to2_controllabel_path'], nifti_template_path, 1, None,  TR)   
    
    else:
        logging.info("Outlier removal: No outliers detected")
        
    logging.info(f"===================================================================")

    return subject