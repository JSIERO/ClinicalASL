"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Main pipeline script for processing ASL MRI data.
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Outlier removal by rejecting dynamics from the data,
    outlier dynamic is one whose mean deltaM (in mask) deviates more than mean deltaM over dynamics + outlier_factor x temporal std brain deltaM 
    precisely the median across the dynamics is used, and mad (median absolute deviation)
    
    Remove outlier dynamic from ASL data based on mean deltaM in GM.
    Parameters:
        deltaM4D: 4D numpy array (x, y, z, time)
        outlier_factor: threshold factor (default 2.5)
        brainmask: brain mask, boolean
        usermask: additional mask, boolean, e.g. gray matter
    Returns:
        deltaM4D with outliers removed: 4D numpy array with outlier dynamics removed
        index_outliers: indices of removed dynamics

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


def save_outlier_removed(data, path_key, context_data, label):
    new_path = append_or(context_data[path_key])
    logging.info(f"Saving ASL outlier-removed data: {label}")
    save_data_nifti(data, new_path, context_data['templateNII_path'], 1, None, context_data['TR'])
    context_data[path_key] = new_path


def asl_outlier_removal(subject, context_tag, usermask=None):
    #
    # Remove outlier dynamics from ASL deltaM signal based on mean deltaM in a combined mask.
    #
    context_data = subject[context_tag]
    outlier_factor = subject['outlier_factor']
    NREPEATS = context_data['NREPEATS']
    NPLDS = context_data['NPLDS']

    # Exclude M0 (first dynamic): shape becomes (x, y, z, NREPEATS, NPLDS, 2)
    ASL_allPLD = context_data['M0ASL_allPLD'][:, :, :, 1:, :, :]

    # Compute deltaM = label - control across time
    deltaM5D = ASL_allPLD[..., 1] - ASL_allPLD[..., 0]  # shape: (x, y, z, NREPEATS, NPLDS)
    deltaM4D = deltaM5D.sum(axis=4)  # shape: (x, y, z, NREPEATS)

    # Load mask and combine with user-supplied mask if provided
    brainmask = nib.load(context_data['mask_path']).get_fdata() > 0
    mask = brainmask if usermask is None else np.logical_and(brainmask, usermask)

    # Masked view of data: shape (N_voxels, NREPEATS)
    voxels_in_mask = deltaM4D[mask]
    voxels_in_mask = np.where(np.isinf(voxels_in_mask), 0, voxels_in_mask)

    # Mean deltaM per dynamic
    mean_per_dyn = np.nanmean(voxels_in_mask, axis=0)

    if np.isnan(mean_per_dyn).any():
        raise ValueError(f"NaNs found in deltaM means at volumes: {np.where(np.isnan(mean_per_dyn))[0]}")

    # Robust outlier detection using MAD
    mu = np.median(mean_per_dyn)
    std = 1.4826 * mad(mean_per_dyn)
    keep_dynamics = np.abs(mean_per_dyn - mu) < outlier_factor * std

    outlier_indices = np.where(~keep_dynamics)[0]
    context_data['outlier_dynamics_indices'] = outlier_indices

    if outlier_indices.size > 0:
        logging.info(f"Outlier removal: Volumes removed (1-based): {(outlier_indices + 1).tolist()}")
    else:
        logging.info("Outlier removal: No outliers detected")

    # Load original PLD data and reshape
    PLDall = nib.load(context_data['PLDall_controllabel_path']).get_fdata()
    PLD_shape = PLDall.shape
    PLDall_reshaped = PLDall.reshape(*PLD_shape[:3], NREPEATS, NPLDS, 2)

    # Keep only the non-outlier repeats
    keep_indices = np.where(keep_dynamics)[0]
    PLDall_kept = PLDall_reshaped[:, :, :, keep_indices, :, :]

    # Reshape back to 4D: (x, y, z, 2 * len(kept) * NPLDS)
    PLDall_or = PLDall_kept.reshape(*PLD_shape[:3], 2 * len(keep_indices) * NPLDS)

    # Derive subsets
    PLD2tolast_or = PLDall_or[:, :, :, NREPEATS*2:]       # for CBF
    PLD1to2_or = PLDall_or[:, :, :, :NREPEATS*2*2]        # for AAT

    # Save all outputs
    save_outlier_removed(PLDall_or, 'PLDall_controllabel_path', context_data, 'all PLDs for AAT')
    save_outlier_removed(PLD2tolast_or, 'PLD2tolast_controllabel_path', context_data, '2-to-last PLDs for CBF')
    save_outlier_removed(PLD1to2_or, 'PLD1to2_controllabel_path', context_data, '1-to-2 PLDs for AAT')

    return subject