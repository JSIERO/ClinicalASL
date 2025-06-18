"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Main pipeline script for processing ASL MRI data.
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Outlier removal by rejecting dynamics/timepoints from the data, based on 
    outlier dynamic is one that has a mean brain CBF (deltaM) > timepoint mean brain CBF (deltaM) + 2.5 temporal std  brain CBF (deltaM)
    precisely the median across the timepoints is used, and mad (median absolute deviation)
    
    Remove outlier timepoints from ASL data based on mean CBF in GM.
    Parameters:
        CBFdata4D: 4D numpy array (x, y, z, time)
        outlierFactor: threshold factor (default 2.5)
        brainmask: brain mask, boolean
        usermask: additional mask, boolean, e.g. gray matter
    Returns:
        CBFdata4D with outliers removed: 4D numpy array with outlier timepoints removed
        index_outliers: indices of removed timepoints

License: BSD 3-Clause License
"""

import numpy as np

def mad(data, axis=None):
    # Median Absolute Deviation: a robust version of standard deviation."""
    return np.median(np.abs(data - np.median(data, axis)), axis)

def asl_outlier_removal(CBFdata4D, brainmask, usermask=None, outlierFactor=2.5):
    #
    # Remove outlier timepoints from ASL data based on mean CBF in a combined mask.
    #
    mask = brainmask if usermask is None else np.logical_and(brainmask, usermask)

    # Create masked 2D view: (voxels, time)
    voxels_in_mask = CBFdata4D[mask]  # shape: (N_voxels, N_timepoints)

    # Replace inf with 0
    voxels_in_mask = np.where(np.isinf(voxels_in_mask), 0, voxels_in_mask)

    # Compute mean CBF per timepoint, ignoring NaNs
    mean_cbf_mask = np.nanmean(voxels_in_mask, axis=0)

    # Check for NaNs
    if np.isnan(mean_cbf_mask).any():
        raise ValueError(f"NaN value found in meanCBF for volumes: {np.where(np.isnan(mean_cbf_mask))[0]}")

    mu_cbf = np.median(mean_cbf_mask)
    std_cbf = 1.4826 * mad(mean_cbf_mask) # median absolute deviation

    keep_mask = np.abs(mean_cbf_mask - mu_cbf) <= outlierFactor * std_cbf
    index_outliers = np.where(~keep_mask)[0]

    if index_outliers.size > 0:
        print(f"Step1: Volume(s) removed (|meanCBF - mu| > {outlierFactor}*std): {index_outliers.tolist()}")
    else:
        print(f"Step1: No outliers detected")

    print("")

    return CBFdata4D[..., keep_mask], index_outliers
