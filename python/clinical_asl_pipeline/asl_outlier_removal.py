"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Main pipeline script for processing ASL MRI data.
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Outlier removal by rejecting dynamics/timepoints from the data, based on 
    outlier dynamic is one that has a mean brain deltaM (deltaM) > timepoint mean brain deltaM (deltaM) + 2.5 temporal std  brain deltaM (deltaM)
    precisely the median across the timepoints is used, and mad (median absolute deviation)
    
    Remove outlier timepoints from ASL data based on mean deltaM in GM.
    Parameters:
        deltaMdata4D: 4D numpy array (x, y, z, time)
        outlier_factor: threshold factor (default 2.5)
        brainmask: brain mask, boolean
        usermask: additional mask, boolean, e.g. gray matter
    Returns:
        deltaMdata4D with outliers removed: 4D numpy array with outlier timepoints removed
        index_outliers: indices of removed timepoints

License: BSD 3-Clause License
"""

import numpy as np
import nibabel as nib

def mad(data, axis=None):
    # Median Absolute Deviation: a robust version of standard deviation."""
    return np.median(np.abs(data - np.median(data, axis)), axis)

def asl_outlier_removal(subject, context_tag, usermask=None):
    #
    # Remove outlier timepoints from ASL deltaM signal based on mean deltaM in a combined mask.
    #
    context_data = subject[context_tag]
    outlier_factor = subject['outlier_factor']

    # M0ASL_allPLD is 6D numpy array (x, y, z,  NDYNS, NPLDS, control/label)
    # ASL_allPLD is 6D numpy array (x, y, z,  NREPEATS, NPLDS, control/label)
    ASL_allPLD = context_data['M0ASL_allPLD'][:,:,:,1:,:,:]  # exclude M0 (first dynamic)

    dims = ASL_allPLD.shape

    deltaM4D = np.zeros((dims[0], dims[1], dims[2], context_data['NREPEATS']))
    deltaM5D = np.zeros((dims[0], dims[1], dims[2], context_data['NREPEATS'],context_data['NPLDS']))

    for i in range(context_data['NREPEATS']):
        for j in range(context_data['NPLDS']):
        # Sum over the PLD dimension, which is axis 4 after selecting the i-th dynamic
            #deltaM4D[:,:,:,i] = np.sum(ASL_allPLD[:,:,:,i,:,1] - ASL_allPLD[:,:,:,i,:,0], axis=3)  # axis=3 is PLD after 
            deltaM5D[:,:,:,i,j] = ASL_allPLD[:,:,:,i,j,1] - ASL_allPLD[:,:,:,i,j,0]
    
    deltaM4D = np.sum(deltaM5D, axis=4)
    import matplotlib.pyplot as plt

    # Plot a mosaic view: middle slice for all 23 dynamics (timepoints)
    mid_slice = deltaM4D.shape[2] // 2
    n_timepoints = deltaM4D.shape[3]
    n_cols = 6
    n_rows = int(np.ceil(n_timepoints / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2*n_cols, 2*n_rows))
    axes = axes.flatten()

    for i in range(n_timepoints):
        axes[i].imshow(deltaM4D[:, :, mid_slice, i], cmap='gray')
        axes[i].set_title(f'Dyn {i+1}')
        axes[i].axis('off')

    # Hide unused subplots
    for j in range(n_timepoints, len(axes)):
        axes[j].axis('off')

    plt.suptitle('Middle Slice Mosaic of deltaM4D (All Dynamics)')
    plt.tight_layout()
    plt.show()
    # Load brainmask
    brainmask = nib.load(context_data['mask_path']).get_fdata()
    brainmask = brainmask > 0  # Ensure binary mask
    
    # combine brainmask with user supplied mask if exists
    mask = brainmask if usermask is None else np.logical_and(brainmask, usermask)

    # Create masked 2D view: (voxels, time)
    voxels_in_mask = deltaM4D[mask]  # shape: (N_voxels, N_timepoints)

    # Replace inf with 0
    voxels_in_mask = np.where(np.isinf(voxels_in_mask), 0, voxels_in_mask)

    # Compute mean deltaM in mask per timepoint, ignoring NaNs
    mean_deltaM_mask = np.nanmean(voxels_in_mask, axis=0)

    # Check for NaNs
    if np.isnan(mean_deltaM_mask).any():
        raise ValueError(f"NaN value found in meandeltaM for volumes: {np.where(np.isnan(mean_deltaM_mask))[0]}")

    mu_deltaM = np.median(mean_deltaM_mask)
    std_deltaM = 1.4826 * mad(mean_deltaM_mask) # median absolute deviation

    keep_mask = np.abs(mean_deltaM_mask - mu_deltaM) <= outlier_factor * std_deltaM
    index_outliers = np.where(~keep_mask)[0]

    if index_outliers.size > 0:
        print(f"Outlier removal: Volume(s) removed (|meandeltaM - mu| > {outlier_factor}*std): {index_outliers.tolist()}")
    else:
        print(f"Outlier removal: No outliers detected")

    return deltaM4D[..., keep_mask], index_outliers
