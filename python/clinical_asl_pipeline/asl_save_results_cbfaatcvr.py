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
import logging
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

    # Load CBF, AAT, ATA from QASL analysis and masks, and registered post- to preACZ NIFTIs
    subject['preACZ']['CBF'] = nib.load(subject['preACZ']['QASL_CBF_path']).get_fdata()
    subject['postACZ']['CBF'] = nib.load(subject['postACZ']['QASL_CBF_path']).get_fdata()
    subject['postACZ']['CBF_2preACZ'] = nib.load(subject['postACZ']['CBF_2preACZ_path']).get_fdata()
    subject['preACZ']['AAT'] = nib.load(subject['preACZ']['QASL_AAT_path']).get_fdata()
    subject['postACZ']['AAT'] = nib.load(subject['postACZ']['QASL_AAT_path']).get_fdata()
    subject['postACZ']['AAT_2preACZ'] = nib.load(subject['postACZ']['AAT_2preACZ_path']).get_fdata()
    subject['preACZ']['ATA'] = nib.load(subject['preACZ']['QASL_ATA_path']).get_fdata()
    subject['postACZ']['ATA'] = nib.load(subject['postACZ']['QASL_ATA_path']).get_fdata()
    subject['postACZ']['ATA_2preACZ'] = nib.load(subject['postACZ']['ATA_2preACZ_path']).get_fdata()
    subject['postACZ']['mask_2preACZ'] = nib.load(subject['postACZ']['mask_2preACZ_path']).get_fdata()
        
    # Create combined masks
    subject['postACZ']['nanmask_2preACZ'] = np.where(subject['postACZ']['mask_2preACZ'], 1.0, np.nan)
    subject['nanmask_combined'] = subject['preACZ']['nanmask'] * subject['postACZ']['nanmask_2preACZ']

    # Compute CVR
    subject['CVR'] = subject['postACZ']['CBF_2preACZ'] - subject['preACZ']['CBF']

    # Smoothing
    subject['CVR_smth'] = asl_smooth_image(subject['CVR'] * subject['nanmask_combined'], 2, subject['FWHM'], subject['preACZ']['VOXELSIZE'])
    subject['preACZ']['AAT_smth'] = asl_smooth_image(subject['preACZ']['AAT'] * subject['preACZ']['nanmask'], 2, subject['FWHM'], subject['postACZ']['VOXELSIZE'])
    subject['postACZ']['AAT_smth'] = asl_smooth_image(subject['postACZ']['AAT'] * subject['postACZ']['nanmask'], 2, subject['FWHM'], subject['postACZ']['VOXELSIZE'])
    subject['postACZ']['AAT_2preACZ_smth'] = asl_smooth_image(subject['postACZ']['AAT_2preACZ'] * subject['postACZ']['nanmask_2preACZ'], 2, subject['FWHM'], subject['postACZ']['VOXELSIZE'])

    # Save results as NIfTI in ASLdir
    save_data_nifti(subject['preACZ']['CBF'], subject['preACZ']['output_CBF_path'], subject['preACZ']['templateNII_path'], 1, None, subject['preACZ']['TR'])
    save_data_nifti(subject['postACZ']['CBF'], subject['postACZ']['output_CBF_path'], subject['postACZ']['templateNII_path'], 1, None, subject['postACZ']['TR'])
    save_data_nifti(subject['postACZ']['CBF_2preACZ'], subject['postACZ']['CBF_2preACZ_path'], subject['preACZ']['templateNII_path'], 1, None, subject['postACZ']['TR'])
    save_data_nifti(subject['preACZ']['AAT'], subject['preACZ']['output_AAT_path'], subject['preACZ']['templateNII_path'], 1, None, subject['preACZ']['TR'])
    save_data_nifti(subject['postACZ']['AAT'], subject['preACZ']['output_AAT_path'], subject['postACZ']['templateNII_path'], 1, None, subject['postACZ']['TR'])
    save_data_nifti(subject['postACZ']['AAT_2preACZ'], subject['postACZ']['AAT_2preACZ_path'], subject['postACZ']['templateNII_path'], 1, None, subject['postACZ']['TR'])
    save_data_nifti(subject['preACZ']['ATA'], subject['preACZ']['output_ATA_path'], subject['preACZ']['templateNII_path'], 1, None, subject['preACZ']['TR'])
    save_data_nifti(subject['postACZ']['ATA'], subject['preACZ']['output_ATA_path'], subject['postACZ']['templateNII_path'], 1, None, subject['postACZ']['TR'])
    save_data_nifti(subject['postACZ']['ATA_2preACZ'], subject['postACZ']['ATA_2preACZ_path'], subject['preACZ']['templateNII_path'], 1, None, subject['postACZ']['TR'])
    save_data_nifti(subject['CVR_smth'], subject['output_CVR_path'], subject['preACZ']['templateNII_path'], 1, None, subject['preACZ']['TR'])

    # Save PNGs
    save_figure_to_png(subject['preACZ']['CBF'], subject['nanmask_combined'], subject['range_cbf'], subject['RESULTSdir'], 'preACZ_CBF', 'CBF', 'viridis')
    save_figure_to_png(subject['postACZ']['CBF_2preACZ'], subject['nanmask_combined'], subject['range_cbf'], subject['RESULTSdir'], 'postACZ_CBF_2preACZ', 'CBF', 'viridis')
    save_figure_to_png(subject['CVR_smth'], subject['nanmask_combined'], subject['range_cvr'], subject['RESULTSdir'], 'CVR', 'CVR', 'vik')
    save_figure_to_png(subject['preACZ']['AAT_smth'], subject['nanmask_combined'], subject['range_AAT'], subject['RESULTSdir'], 'preACZ_AAT', 'time', 'devon')
    save_figure_to_png(subject['postACZ']['AAT_2preACZ_smth'], subject['nanmask_combined'], subject['range_AAT'], subject['RESULTSdir'], 'postACZ_AAT_2preACZ', 'time', 'devon')
    save_figure_to_png(subject['preACZ']['ATA'], subject['nanmask_combined'], subject['range_ATA'], subject['RESULTSdir'], 'preACZ_ATA', 'CBF', 'viridis')
    save_figure_to_png(subject['postACZ']['ATA_2preACZ'], subject['nanmask_combined'], subject['range_ATA'], subject['RESULTSdir'], 'postACZ_ATA_2preACZ', 'CBF', 'viridis')

    save_data_dicom(subject['preACZ']['CBF'], os.path.join(subject['DICOMsubjectdir'], subject['preACZ']['templateDCM_CBF_path']), os.path.join(subject['DICOMoutputdir'], subject['preACZ']['templateDCM_CBF_path'] + '.dcm'), 'WIP preACZ CBF MD-ASL', subject['range_cbf'], 'CBF')
    save_data_dicom(subject['preACZ']['AAT_smth'], os.path.join(subject['DICOMsubjectdir'], subject['preACZ']['templateDCM_AAT_path']), os.path.join(subject['DICOMoutputdir'], subject['preACZ']['templateDCM_AAT_path'] + '.dcm'), 'WIP preACZ AAT(s) MD-ASL', subject['range_AAT'], 'AAT')
    save_data_dicom(subject['preACZ']['ATA'], os.path.join(subject['DICOMsubjectdir'], subject['preACZ']['templateDCM_ATA_path']), os.path.join(subject['DICOMoutputdir'], subject['preACZ']['templateDCM_ATA_path'] + '.dcm'), 'WIP preACZ ATA MD-ASL', subject['range_ATA'], 'ATA')
    save_data_dicom(subject['CVR_smth'], os.path.join(subject['DICOMsubjectdir'], subject['preACZ']['templateDCM_CVR_path']), os.path.join(subject['DICOMoutputdir'], subject['preACZ']['templateDCM_CVR_path'] + '.dcm'), 'WIP CVR MD-ASL', subject['range_cvr'], 'CVR')
    
    save_data_dicom(subject['postACZ']['CBF'], os.path.join(subject['DICOMsubjectdir'], subject['postACZ']['templateDCM_CBF_path']), os.path.join(subject['DICOMoutputdir'], subject['postACZ']['templateDCM_CBF_path'] + '.dcm'), 'WIP postACZ CBF MD-ASL', subject['range_cbf'], 'CBF')
    save_data_dicom(subject['postACZ']['AAT_smth'], os.path.join(subject['DICOMsubjectdir'], subject['postACZ']['templateDCM_AAT_path']), os.path.join(subject['DICOMoutputdir'], subject['postACZ']['templateDCM_AAT_path'] + '.dcm'), 'WIP postACZ AAT(s) MD-ASL', subject['range_AAT'], 'AAT')
    save_data_dicom(subject['postACZ']['ATA'], os.path.join(subject['DICOMsubjectdir'], subject['postACZ']['templateDCM_ATA_path']), os.path.join(subject['DICOMoutputdir'], subject['postACZ']['templateDCM_ATA_path'] + '.dcm'), 'WIP postACZ ATA MD-ASL', subject['range_ATA'], 'ATA')

    logging.info(f"Results:  PACS-ready DICOMS, NIFTI, .png's saved for subject {subject['ASLdir']}")
    
