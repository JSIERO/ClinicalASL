"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    This script saves CBF, AAT, ATA, and CVR results for a subject after ASL processing.
    It loads registered and original NIfTI images, creates combined masks, computes CVR,
    applies smoothing, and saves results as NIfTI, PNG, and DICOM files for further analysis and PACS export.

License: BSD 3-Clause License
"""

import os
import numpy as np
import nibabel as nib
from clinical_asl_pipeline.asl_smooth_image import asl_smooth_image
from clinical_asl_pipeline.utils.save_figure_to_png import save_figure_to_png  
from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti
from clinical_asl_pipeline.utils.save_data_dicom import save_data_dicom

def asl_save_results_cbfaatcvr(subject):
# Save CBF, AAT, ATA, and CVR results for a subject.
# Steps:
# 1. Load registered and original NIfTI images into subject dict.
# 2. Create combined masks for valid voxels.
# 3. Compute CVR as difference between postACZ and preACZ CBF.
# 4. Smooth CVR, AAT images using Gaussian smoothing.
# 5. Save results as NIfTI files.
# 6. Save visualizations as PNG images.
# 7. Save results as DICOM files for Philips multi-frame export for PACS

    # Load CBF, AAT, masks, and registered post- to preACZ NIFTIs
    subject['preACZ']['CBF'] = nib.load(subject['preACZ_CBF_path']).get_fdata()
    subject['postACZ']['CBF'] = nib.load(subject['postACZ_CBF_path']).get_fdata()
    subject['postACZ']['CBF_2preACZ'] = nib.load(subject['postACZ_CBF_2preACZ_path']).get_fdata()
    subject['preACZ']['AAT'] = nib.load(subject['preACZ_AAT_path']).get_fdata()
    subject['postACZ']['AAT'] = nib.load(subject['postACZ_AAT_path']).get_fdata()
    subject['postACZ']['AAT_2preACZ'] = nib.load(subject['postACZ_AAT_2preACZ_path']).get_fdata()
    subject['preACZ']['ATA'] = nib.load(subject['preACZ_ATA_path']).get_fdata()
    subject['postACZ']['ATA'] = nib.load(subject['postACZ_ATA_path']).get_fdata()
    subject['postACZ']['ATA_2preACZ'] = nib.load(subject['postACZ_ATA_2preACZ_path']).get_fdata()
    subject['postACZ']['mask_2preACZ'] = nib.load(subject['postACZ_mask_2preACZ_path']).get_fdata()
        
    # Create combined masks
    subject['postACZ']['nanmask_2preACZ'] = np.where(subject['postACZ']['mask_2preACZ'], 1.0, np.nan)
    subject['nanmask_combined'] = subject['preACZ']['nanmask'] * subject['postACZ']['nanmask_2preACZ']

    # Compute CVR
    subject['CVR'] = subject['postACZ']['CBF_2preACZ'] - subject['preACZ']['CBF']

    # Smoothing
    subject['CVR_smth'] = asl_smooth_image(subject['CVR'] * subject['nanmask_combined'], 2, subject['FWHM'], subject['VOXELSIZE'])
    subject['preACZ']['AAT_smth'] = asl_smooth_image(subject['preACZ']['AAT'] * subject['preACZ']['nanmask'], 2, subject['FWHM'], subject['VOXELSIZE'])
    subject['postACZ']['AAT_smth'] = asl_smooth_image(subject['postACZ']['AAT'] * subject['postACZ']['nanmask'], 2, subject['FWHM'], subject['VOXELSIZE'])
    subject['postACZ']['AAT_2preACZ_smth'] = asl_smooth_image(subject['postACZ']['AAT_2preACZ'] * subject['postACZ']['nanmask_2preACZ'], 2, subject['FWHM'], subject['VOXELSIZE'])

    # Save results as NIfTI in ASLdir
    save_data_nifti(subject['preACZ']['CBF'], os.path.join(subject['ASLdir'], 'preACZ_CBF.nii.gz'), subject[f'preACZ_templateNII_path'], 1, None, subject['TR'])
    save_data_nifti(subject['postACZ']['CBF'], os.path.join(subject['ASLdir'], 'postACZ_CBF.nii.gz'), subject[f'postACZ_templateNII_path'], 1, None, subject['TR'])
    save_data_nifti(subject['postACZ']['CBF_2preACZ'], os.path.join(subject['ASLdir'], 'postACZ_CBF_2preACZ.nii.gz'), subject[f'preACZ_templateNII_path'], 1, None, subject['TR'])
    save_data_nifti(subject['preACZ']['AAT'], os.path.join(subject['ASLdir'], 'preACZ_AAT.nii.gz'), subject[f'preACZ_templateNII_path'], 1, None, subject['TR'])
    save_data_nifti(subject['postACZ']['AAT'], os.path.join(subject['ASLdir'], 'postACZ_AAT.nii.gz'), subject[f'postACZ_templateNII_path'], 1, None, subject['TR'])
    save_data_nifti(subject['postACZ']['AAT_2preACZ'], os.path.join(subject['ASLdir'], 'postACZ_AAT_2preACZ.nii.gz'), subject[f'postACZ_tpreACZ_templateNII_pathemplateNII_path'], 1, None, subject['TR'])
    save_data_nifti(subject['preACZ']['ATA'], os.path.join(subject['ASLdir'], 'preACZ_ATA.nii.gz'), subject[f'preACZ_templateNII_path'], 1, None, subject['TR'])
    save_data_nifti(subject['postACZ']['ATA'], os.path.join(subject['ASLdir'], 'postACZ_ATA.nii.gz'), subject[f'postACZ_templateNII_path'], 1, None, subject['TR'])
    save_data_nifti(subject['postACZ']['ATA_2preACZ'], os.path.join(subject['ASLdir'], 'postACZ_ATA_2preACZ.nii.gz'), subject[f'preACZ_templateNII_path'], 1, None, subject['TR'])
    save_data_nifti(subject['CVR_smth'], os.path.join(subject['ASLdir'], 'CVR_smth.nii.gz'), subject[f'preACZ_templateNII_path'], 1, None, subject['TR'])

    # Save PNGs
    save_figure_to_png(subject['preACZ']['CBF'], subject['nanmask_combined'], subject['range_cbf'], subject['RESULTSdir'], f"preACZ_CBF", 'CBF', 'viridis')
    save_figure_to_png(subject['postACZ']['CBF_2preACZ'], subject['nanmask_combined'], subject['range_cbf'], subject['RESULTSdir'], f"postACZ_CBF_2preACZ", 'CBF', 'viridis')
    save_figure_to_png(subject['CVR_smth'], subject['nanmask_combined'], subject['range_cvr'], subject['RESULTSdir'], 'CVR_smth', 'CVR', 'vik')
    save_figure_to_png(subject['preACZ']['AAT_smth'], subject['nanmask_combined'], subject['range_AAT'], subject['RESULTSdir'], 'preACZ_AAT_smth', 'time', 'devon')
    save_figure_to_png(subject['postACZ']['AAT_2preACZ_smth'], subject['nanmask_combined'], subject['range_AAT'], subject['RESULTSdir'], 'postACZ_AAT_2preACZ_smth', 'time', 'devon')
    save_figure_to_png(subject['preACZ']['ATA'], subject['nanmask_combined'], subject['range_ATA'], subject['RESULTSdir'], 'preACZ_ATA', 'CBF', 'viridis')
    save_figure_to_png(subject['postACZ']['ATA_2preACZ'], subject['nanmask_combined'], subject['range_ATA'], subject['RESULTSdir'], 'postACZ_ATA_2preACZ', 'CBF', 'viridis')

    print("Saving DICOM files...")
    save_data_dicom(subject['preACZ']['CBF'], os.path.join(subject['DICOMdir'], subject['preACZfilenameDCM_CBF']), os.path.join(subject['DICOMRESULTSdir'], subject['preACZfilenameDCM_CBF'] + '.dcm'), 'WIP preACZ CBF MD-ASL', subject['range_cbf'], 'CBF')
    save_data_dicom(subject['preACZ']['AAT_smth'], os.path.join(subject['DICOMdir'], subject['preACZfilenameDCM_AAT']), os.path.join(subject['DICOMRESULTSdir'], subject['preACZfilenameDCM_AAT'] + '.dcm'), 'WIP preACZ AAT(s) MD-ASL', subject['range_AAT'], 'AAT')
    save_data_dicom(subject['preACZ']['ATA'], os.path.join(subject['DICOMdir'], subject['preACZfilenameDCM_ATA']), os.path.join(subject['DICOMRESULTSdir'], subject['preACZfilenameDCM_ATA'] + '.dcm'), 'WIP preACZ ATA MD-ASL', subject['range_ATA'], 'ATA')
    
    save_data_dicom(subject['CVR_smth'], os.path.join(subject['DICOMdir'], subject['preACZfilenameDCM_CVR']), os.path.join(subject['DICOMRESULTSdir'], subject['preACZfilenameDCM_CVR'] + '.dcm'), 'WIP CVR MD-ASL', subject['range_cvr'], 'CVR')
    save_data_dicom(subject['postACZ']['CBF'], os.path.join(subject['DICOMdir'], subject['postACZfilenameDCM_CBF']), os.path.join(subject['DICOMRESULTSdir'], subject['postACZfilenameDCM_CBF'] + '.dcm'), 'WIP postACZ CBF MD-ASL', subject['range_cbf'], 'CBF')
    save_data_dicom(subject['postACZ']['AAT_smth'], os.path.join(subject['DICOMdir'], subject['postACZfilenameDCM_AAT']), os.path.join(subject['DICOMRESULTSdir'], subject['postACZfilenameDCM_AAT'] + '.dcm'), 'WIP postACZ AAT(s) MD-ASL', subject['range_AAT'], 'AAT')
    save_data_dicom(subject['postACZ']['ATA'], os.path.join(subject['DICOMdir'], subject['postACZfilenameDCM_ATA']), os.path.join(subject['DICOMRESULTSdir'], subject['postACZfilenameDCM_ATA'] + '.dcm'), 'WIP postACZ ATA MD-ASL', subject['range_ATA'], 'ATA')

    print("CBF, AAT, ATA, CVR Results: DICOMs and PNGs created")

