"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL T1-from-M0 computation module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Functions for brain extraction, computing T1-weighted images from multi-PLD M0 data, and tissue segmentation.

License: BSD 3-Clause License
"""

import logging
import numpy as np
import nibabel as nib

from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti

def asl_t1_from_m0(subject, context_tag):  
    
    # Process M0 image to compute T1w image from multi-PLD M0 data.
    #
    # This includes brain extraction, Look-Locker correction, and T1w from mPLD M0 computation.
    #
    # subject:  Dictionary containing subject data and parameters.
    # context_tag:      String context_tag for the subject's data, eg 'baseline', 'stimulus', etc.
    #
    # Returns the updated subject dictionary with T1w image and masks.
    # Ensure that the subject dictionary contains necessary keys like 'ASLdir', 'PLDS', 'NPLDS', 'LookLocker_correction_factor_perPLD', etc.
    # Example usage:      
    # subject = asl_t1_from_m0(subject, 'baseline')  for baseline ASL data before diamox #  
    #       

    # Use a shorter alias for context_data
    context_data = subject[context_tag]     
    
    LookLocker_correction_factor_perPLD = context_data['LookLocker_correction_factor_perPLD']
    T1fromM0_path = context_data['T1fromM0_path']
    templateNIFTI_path =context_data['templateNIFTI_path']
    M0_allPLD = context_data['M0_allPLD']
    brainmask = context_data['mask'] 
    PLDS =  context_data['PLDS']
    TR = context_data['TR']

    # Remove the Look-Locker correction (by multiplication) to compute the T1w profile (vectorized)
    M0_allPLD_noLLcorr = M0_allPLD * LookLocker_correction_factor_perPLD[None, None, None, :]
    logging.info('ASL T1-from-M0 computation started')

    # Compute T1w image from multi-PLD M0 and save nifti
    logging.info('Create T1w image from multiPLD M0')   

    T1fromM0 = asl_t1_from_m0_compute(M0_allPLD_noLLcorr, brainmask, PLDS)
    save_data_nifti(T1fromM0, T1fromM0_path, templateNIFTI_path, 1, [0, 500], TR)

    # Final T1fromM0 load
    context_data['T1fromM0'] = nib.load(T1fromM0_path).get_fdata()
    return subject

def asl_t1_from_m0_compute(DATA4D, MASK, TIMEARRAY):
    # Compute T1w image from M0 data using a linear fit.  "  
    # DATA4D: 4D numpy array of M0 data (x, y, z, PLD)
    # MASK: 3D numpy array of brain mask (x, y, z)
    # TIMEARRAY: 1D numpy array of PLD times in seconds
    DATA4D_shape = DATA4D.shape
    DATA2D = DATA4D.reshape(-1, DATA4D_shape[3])
    brain_voxels = MASK.flatten() > 0

    with np.errstate(divide='ignore'):
        DATA2D_log = np.log(DATA2D)

    CONSTANTARRAY = np.ones((DATA2D_log.shape[1], 1))
    T1fit = np.zeros((DATA2D_log.shape[0], 2))
    b = TIMEARRAY.reshape(-1, 1)

    for i in range(DATA2D_log.shape[0]):
        if not brain_voxels[i]:
            continue

        A_row = DATA2D_log[i, :]
        if np.any(np.isnan(A_row)) or np.any(np.isinf(A_row)):
            continue

        A = np.hstack([CONSTANTARRAY, A_row.reshape(-1, 1)])
        try:
            fit_result = np.linalg.lstsq(A, b, rcond=None)[0]
            T1fit[i, :] = fit_result.flatten()
        except np.linalg.LinAlgError:
            T1fit[i, :] = [0, 0]

    data_R1fit = T1fit[:, 1].reshape(DATA4D_shape[0], DATA4D_shape[1], DATA4D_shape[2])

    with np.errstate(divide='ignore', invalid='ignore'):
        data_T1fit_brain = (-1 / data_R1fit) * MASK * 1e3
        data_T1fit_brain[~np.isfinite(data_T1fit_brain)] = 0  # Clean NaNs and Infs
        data_T1fit_brain = abs(data_T1fit_brain) # take magnitude of the data so to avoid negative values
    
    data_T1fit_brain_1D = data_T1fit_brain.flatten()[brain_voxels]
    THRESHOLDMASK_FACTOR = np.median(data_T1fit_brain_1D) + 2 * np.std(data_T1fit_brain_1D) # take median + 2 x std as threshold to remove extremely high values
    
    valid_range_mask = (data_T1fit_brain > 0) & (data_T1fit_brain <= THRESHOLDMASK_FACTOR)
    T1fromM0 = data_T1fit_brain * valid_range_mask
    T1fromM0[np.isnan(T1fromM0)] = 0

    return T1fromM0


