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

import os
import numpy as np
import nibabel as nib
import subprocess
from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti

def asl_bet_t1_from_m0(subject, prefix, fast):  
    
    # Process M0 image to compute T1w image from multi-PLD M0 data.
    #
    # This includes brain extraction, Look-Locker correction, and T1w from mPLD M0 computation.
    #
    # subject:  Dictionary containing subject data and parameters.
    # prefix:   String prefix for the subject's data.
    # fast:     String indicating whether to use FSL's FAST for tissue segmentation.
    #
    # Returns the updated subject dictionary with T1w image and masks.
    # Ensure that the subject dictionary contains necessary keys like 'ASLdir', 'PLDS', 'NPLDS', 'LookLocker_correction_factor_perPLD', etc.
    # Example usage:      
    # subject = asl_t1_from_m0(subject, 'preACZ', 'fast')  for baseline ASL data before diamox #  
    #       

    if prefix not in subject:
        subject[prefix] = {}  # Ensure key exists before assigning
     
    # Brain extraction on M0 image using FSL's BET
    m0 = subject[f'{prefix}_m0_path']
    mask = subject[f'{prefix}_mask_path']
    m0_brain_path = os.path.join(subject['ASLdir'], f'{prefix}_M0_brain')

    print('Brain masking: Running BET on M0 image:', mask)   
    subprocess.run(f'bet {m0} {m0_brain_path} -m -f 0.4 -g 0', shell=True, check=True) #

    # Load and process brain mask, binarize it, and create a NaN mask
    mask_data = nib.load(mask).get_fdata()
    subject[prefix]['mask'] = mask_data > 0
    subject[prefix]['nanmask'] = np.where(subject[prefix]['mask'], 1.0, np.nan)

    # Remove the Look-Locker correction (by multiplication) to compute the T1w profile (vectorized)
    M0_allPLD_noLLcorr = subject[prefix]['M0_allPLD'] * subject['LookLocker_correction_factor_perPLD'][None, None, None, :]

    # Compute T1w image from multi-PLD M0 and save nifti
    print('Create T1w image from multiPLD M0')   

    T1fromM0 = asl_t1_from_m0_compute(M0_allPLD_noLLcorr, subject[prefix]['mask'], subject['PLDS'])
    save_data_nifti(T1fromM0, subject[f'{prefix}_T1fromM0_path'], subject[f'{prefix}_templateNII_path'], 1, [0, 500], subject['TR'])

    if fast != 'fast':
        # Segment tissue using FSL FAST
        t1_path = subject[f'{prefix}_T1fromM0_path']
        subprocess.run(f'fast -b -g -B {t1_path}', shell=True)
        subprocess.run(f'fslmaths {t1_path}_restore {t1_path}', shell=True)

        # Load tissue segmentations
        subject[prefix]['CSFmask'] = nib.load(f'{t1_path}_seg_0.nii.gz').get_fdata()
        subject[prefix]['GMmask'] = nib.load(f'{t1_path}_seg_1.nii.gz').get_fdata()
        subject[prefix]['WMmask'] = nib.load(f'{t1_path}_seg_2.nii.gz').get_fdata()

    # Final T1fromM0 load
    subject[prefix]['T1fromM0'] = nib.load(subject[f'{prefix}_T1fromM0_path']).get_fdata()
    return subject

def asl_t1_from_m0_compute(DATA4D, MASK, TIMEARRAY):
    # Compute T1w image from M0 data using a linear fit.  "  
    # DATA4D: 4D numpy array of M0 data (x, y, z, PLD)
    # MASK: 3D numpy array of brain mask (x, y, z)
    # TIMEARRAY: 1D numpy array of PLD times in seconds
    dims = DATA4D.shape
    DATA2D = DATA4D.reshape(-1, dims[3])
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

    data_R1fit = T1fit[:, 1].reshape(dims[0], dims[1], dims[2])

    with np.errstate(divide='ignore', invalid='ignore'):
        data_T1fit_brain = (-1 / data_R1fit) * MASK * 1e3
        data_T1fit_brain[~np.isfinite(data_T1fit_brain)] = 0  # Clean NaNs and Infs

    valid_range_mask = (data_T1fit_brain > 0) & (data_T1fit_brain <= 300)
    T1fromM0 = data_T1fit_brain * valid_range_mask
    T1fromM0[np.isnan(T1fromM0)] = 0

    return T1fromM0


