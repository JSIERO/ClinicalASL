import os
import numpy as np
import nibabel as nib
import subprocess
from clinical_asl_pipeline.save_data_nifti import  save_data_nifti
from clinical_asl_pipeline.asl_smooth_image import  asl_smooth_image

def asl_t1_from_m0_processing(subject, prefix, fast):       
   
    # Brain extraction on M0 image using FSL's BET   
    m0_path = os.path.join(subject['ASLdir'], f"{prefix}_M0")
    brain_path = os.path.join(subject['ASLdir'], f"{prefix}_M0_brain")
    mask_path = os.path.join(subject['ASLdir'], f"{prefix}_M0_brain_mask")
    print(mask_path)
    
    if prefix not in subject:
        subject[prefix] = {}  # Ensure key exists before assigning
    
    subprocess.run(f"bet {m0_path} {brain_path} -m -f 0.4 -g 0", shell=True, check=True)
    output1 = os.path.join(subject['ASLdir'], f"{prefix}_temp_M0_brain_mask_dilM")
    subprocess.run(f"fslmaths {mask_path} -dilM {output1}", shell=True, check=True)
    output2 = os.path.join(subject['ASLdir'], f"{prefix}_temp_M0_brain_mask_ero")
    subprocess.run(f"fslmaths {mask_path} -kernel 2D -ero {output2}", shell=True, check=True)
    output3 = os.path.join(subject['ASLdir'], f"{prefix}_temp_M0_dilD")
    subprocess.run(f"fslmaths {m0_path} -mul {output2} -kernel 2D -dilD {output3}", shell=True, check=True)

    # Load and process brain mask
    brainmask_nii = nib.load(f"{mask_path}.nii.gz")
    brainmask_data = brainmask_nii.get_fdata().astype(bool)
    subject[prefix]['brainmask'] = brainmask_data

    nanmask = brainmask_data.astype(float)
    nanmask[nanmask == 0] = np.nan
    subject[prefix]['nanmask'] = nanmask

    # Load dilated M0 image
    M0_dilD_path = os.path.join(subject['ASLdir'], f"{prefix}_temp_M0_dilD.nii.gz")
    M0_dilD = nib.load(M0_dilD_path).get_fdata()

    if fast != 'fast':
        # Apply smoothing 
        subject[prefix]['M0_forQCBF'] = asl_smooth_image(M0_dilD, 3, subject['FWHM_M0'], subject['VOXELSIZE'])

        # Save smoothed M0 image
        save_data_nifti(subject[prefix]['M0_forQCBF'], os.path.join(subject['ASLdir'], f"{prefix}_M0_forQCBF.nii.gz"), subject['dummyfilenameSaveNII'], 1, [0, 500], subject['TR'])

    # Clean up temp files
    subprocess.run(f"rm {subject['ASLdir']}/*temp*", shell=True)

    # Remove Look-Locker correction for all PLDs
    M0_allPLD_noLLcorr = []
   
    for i in range(subject['NPLDS']):
        corrected = subject[prefix]['M0_allPLD'][:, :, :, i] * subject['LookLocker_correction_factor_perPLD'][i]
        M0_allPLD_noLLcorr.append(corrected)
    M0_allPLD_noLLcorr = np.stack(M0_allPLD_noLLcorr, axis=-1)

    # Compute T1w image from multi-PLD M0
    T1fromM0 = asl_t1_from_m0_compute(M0_allPLD_noLLcorr, subject[prefix]['brainmask'], subject['PLDS'])

    save_data_nifti(T1fromM0, os.path.join(subject['ASLdir'], f"{prefix}_T1fromM0.nii.gz"), subject['dummyfilenameSaveNII'], 1, [0, 500], subject['TR'])

    if fast != 'fast':
        # Segment tissue using FSL FAST
        t1_path = os.path.join(subject['ASLdir'], f"{prefix}_T1fromM0")
        subprocess.run(f"fast -b -g -B {t1_path}", shell=True)
        subprocess.run(f"fslmaths {t1_path}_restore {t1_path}", shell=True)

        # Load tissue segmentations
        subject[prefix]['CSFmask'] = nib.load(f"{t1_path}_seg_0.nii.gz").get_fdata()
        subject[prefix]['GMmask'] = nib.load(f"{t1_path}_seg_1.nii.gz").get_fdata()
        subject[prefix]['WMmask'] = nib.load(f"{t1_path}_seg_2.nii.gz").get_fdata()

    # Final T1fromM0 load
    subject[prefix]['T1fromM0'] = nib.load(os.path.join(subject['ASLdir'],f"{prefix}_T1fromM0.nii.gz")).get_fdata()
    return subject

def asl_t1_from_m0_compute(DATA4D, MASK, TIMEARRAY):
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


