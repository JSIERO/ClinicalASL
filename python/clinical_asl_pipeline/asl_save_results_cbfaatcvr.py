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
        
    # === Load data ===
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

    # === Mask prep ===
    subject['postACZ']['nanmask_2preACZ'] = np.where(subject['postACZ']['mask_2preACZ'], 1.0, np.nan)
    subject['nanmask_combined'] = subject['preACZ']['nanmask'] * subject['postACZ']['nanmask_2preACZ']

    # === Compute CVR ===
    subject['CVR'] = subject['postACZ']['CBF_2preACZ'] - subject['preACZ']['CBF']

    # === Apply smoothing ===
    subject['CVR_smth'] = asl_smooth_image(subject['CVR'] * subject['nanmask_combined'], 2, subject['FWHM'], subject['preACZ']['VOXELSIZE'])
    subject['preACZ']['AAT_smth'] = asl_smooth_image(subject['preACZ']['AAT'] * subject['preACZ']['nanmask'], 2, subject['FWHM'], subject['postACZ']['VOXELSIZE'])
    subject['postACZ']['AAT_smth'] = asl_smooth_image(subject['postACZ']['AAT'] * subject['postACZ']['nanmask'], 2, subject['FWHM'], subject['postACZ']['VOXELSIZE'])
    subject['postACZ']['AAT_2preACZ_smth'] = asl_smooth_image(subject['postACZ']['AAT_2preACZ'] * subject['postACZ']['nanmask_2preACZ'], 2, subject['FWHM'], subject['postACZ']['VOXELSIZE'])

    # === Define output lists ===
    fields_main = [
        ('CBF', 'range_cbf', 'CBF', 'viridis', 'output_CBF_path'),
        ('AAT_smth', 'range_AAT', 'AAT', 'devon', 'output_AAT_path'),
        ('ATA', 'range_ATA', 'ATA', 'viridis', 'output_ATA_path'),
    ]
    fields_cvr = [
        ('CVR_smth', 'range_cvr', 'CVR', 'vik', 'output_CVR_path'),
    ]
    fields_2preACZ = [
        ('CBF_2preACZ', 'range_cbf', 'CBF', 'viridis', 'CBF_2preACZ_path'),
        ('AAT_2preACZ_smth', 'range_AAT', 'AAT', 'devon', 'AAT_2preACZ_path'),
        ('ATA_2preACZ', 'range_ATA', 'ATA', 'viridis', 'ATA_2preACZ_path'),
    ]

    # === Helper: Save NIfTI + DICOM ===
    def save_nifti_and_dicom(phase, fields, allow_dicom=True):
        for field, range_key, label, _, output_key in fields:
            data = subject[phase].get(field) if field != 'CVR_smth' else subject['CVR_smth']
            path = subject[phase].get(output_key) if field != 'CVR_smth' else subject['output_CVR_path']
            template = subject[phase]['templateNII_path'] if field != 'CVR_smth' else subject['preACZ']['templateNII_path']
            TR = subject[phase]['TR'] if field != 'CVR_smth' else subject['preACZ']['TR']

            if data is not None and path:
                save_data_nifti(data, path, template, 1, None, TR)
                logging.info(f"Saved NIfTI: {path}")

                if allow_dicom:
                    dcm_template = subject[phase].get(f'templateDCM_{label}_path')
                    if dcm_template:
                        dcm_outpath = os.path.join(subject['DICOMoutputdir'], dcm_template + '.dcm')
                        save_data_dicom(data,
                                        os.path.join(subject['DICOMsubjectdir'], dcm_template),
                                        dcm_outpath,
                                        f'WIP {phase} {label} MD-ASL',
                                        subject[range_key], label)
                        logging.info(f"Saved DICOM: {dcm_outpath}")

    # === Helper: Save PNG ===
    def save_png(phase, fields):
        for field, range_key, label, cmap, _ in fields:
            data = subject[phase].get(field) if field != 'CVR_smth' else subject['CVR_smth']
            if data is not None:
                png_name = f'{phase}_{label}'
                save_figure_to_png(data, subject['nanmask_combined'], subject[range_key],
                                   subject['RESULTSdir'], png_name, label, cmap)
                logging.info(f"Saved PNG: {png_name}.png")

    # === Execute saves ===
    save_nifti_and_dicom('preACZ', fields_main)
    save_nifti_and_dicom('preACZ', fields_cvr)
    save_nifti_and_dicom('postACZ', fields_main)
    save_nifti_and_dicom('postACZ', fields_2preACZ, allow_dicom=False)

    save_png('preACZ', fields_main + fields_cvr)
    save_png('postACZ', fields_2preACZ)

    # === Final log ===
    logging.info(f"Results complete: PACS-ready DICOMS, NIFTI, .png's saved for subject {subject['ASLdir']}")
