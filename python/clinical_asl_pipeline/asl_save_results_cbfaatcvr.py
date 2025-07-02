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

    Main steps:
    - Load CBF, AAT, ATA, and mask images for both baseline and stimulus conditions.
    - Prepare nanmasks to handle non-brain voxels and combine baseline/stimulus masks.
    - Compute CVR (Cerebrovascular Reactivity) as the difference between stimulus and baseline CBF.
    - Apply spatial smoothing to CVR and AAT maps, with method-dependent logic for AAT smoothing.
    - Define output fields for saving (CBF, AAT, ATA, CVR) and their visualization parameters.
    - Save results in three formats:
        * NIfTI: for quantitative analysis and further processing.
        * DICOM: for clinical PACS integration, using appropriate tags and templates.
        * PNG: for quick visualization and quality control.
    - Modular helper functions are used for saving in each format, with logging for traceability.

    This script is intended to be called as part of the ClinicalASL pipeline after ASL quantification
    and registration steps are complete.

License: BSD 3-Clause License
"""

import os
import logging
import numpy as np
import nibabel as nib
from clinical_asl_pipeline.asl_smooth_image import asl_smooth_image
from clinical_asl_pipeline.utils.save_figure_to_png import save_figure_to_png
from clinical_asl_pipeline.utils.save_png_to_dicom import save_png_to_dicom  
from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti
from clinical_asl_pipeline.utils.save_data_dicom import save_data_dicom

def asl_save_results_cbfaatcvr(subject):
    # === Load data ===
    subject['baseline']['CBF'] = nib.load(subject['baseline']['QASL_CBF_path']).get_fdata()
    subject['stimulus']['CBF'] = nib.load(subject['stimulus']['QASL_CBF_path']).get_fdata()
    subject['stimulus']['CBF_2baseline'] = nib.load(subject['stimulus']['CBF_2baseline_path']).get_fdata()
    subject['baseline']['AAT'] = nib.load(subject['baseline']['QASL_AAT_path']).get_fdata()
    subject['stimulus']['AAT'] = nib.load(subject['stimulus']['QASL_AAT_path']).get_fdata()
    subject['stimulus']['AAT_2baseline'] = nib.load(subject['stimulus']['AAT_2baseline_path']).get_fdata()
    subject['baseline']['ATA'] = nib.load(subject['baseline']['QASL_ATA_path']).get_fdata()
    subject['stimulus']['ATA'] = nib.load(subject['stimulus']['QASL_ATA_path']).get_fdata()
    subject['stimulus']['ATA_2baseline'] = nib.load(subject['stimulus']['ATA_2baseline_path']).get_fdata()
    subject['stimulus']['mask_2baseline'] = nib.load(subject['stimulus']['mask_2baseline_path']).get_fdata()

    # === Mask prep ===
    subject['stimulus']['nanmask_2baseline'] = np.where(subject['stimulus']['mask_2baseline'], 1.0, np.nan)
    subject['nanmask_combined'] = subject['baseline']['nanmask'] * subject['stimulus']['nanmask_2baseline']

    # === Compute CVR ===
    subject['CVR'] = subject['stimulus']['CBF_2baseline'] - subject['baseline']['CBF']

    # === Apply smoothing, use nanmask to preserve outside brain edges ===
    subject['CVR_smth'] = asl_smooth_image(subject['CVR'] * subject['nanmask_combined'], 2, subject['FWHM'], subject['baseline']['VOXELSIZE'])

    # do not smooth AAT when using QASL SSVB (inherently does spatial smoothing)
    if subject['inference_method'] == 'ssvb':  
        subject['baseline']['AAT_smth'] = subject['baseline']['AAT'] 
        subject['stimulus']['AAT_smth'] = subject['stimulus']['AAT']
        subject['stimulus']['AAT_2baseline_smth'] = subject['stimulus']['AAT_2baseline']
    # do smooth AAT when using QASL VABY or BASIl 
    else:
        subject['baseline']['AAT_smth'] = asl_smooth_image(subject['baseline']['AAT'] * subject['baseline']['nanmask'], 2, subject['FWHM'], subject['stimulus']['VOXELSIZE'])
        subject['stimulus']['AAT_smth'] = asl_smooth_image(subject['stimulus']['AAT'] * subject['stimulus']['nanmask'], 2, subject['FWHM'], subject['stimulus']['VOXELSIZE'])
        subject['stimulus']['AAT_2baseline_smth'] = asl_smooth_image(subject['stimulus']['AAT_2baseline'] * subject['stimulus']['nanmask_2baseline'], 2, subject['FWHM'], subject['stimulus']['VOXELSIZE'])   

    # === Define output lists ===
    # for all the context data, baseline, stimulus


    fields_main = [
        ('CBF', 'range_cbf', 'CBF', 'viridis', 'output_CBF_path'),
        ('AAT_smth', 'range_AAT', 'AAT', 'devon', 'output_AAT_path'),
        ('ATA', 'range_ATA', 'ATA', 'viridis', 'output_ATA_path'),
    ]
    fields_cvr = [
        ('CVR_smth', 'range_cvr', 'CVR', 'vik', 'output_CVR_path'),
    ]
    fields_2baseline = [ # for the registered stimulus output to baseline
        ('CBF_2baseline', 'range_cbf', 'CBF', 'viridis', 'CBF_2baseline_path'),
        ('AAT_2baseline_smth', 'range_AAT', 'AAT', 'devon', 'AAT_2baseline_path'),
        ('ATA_2baseline', 'range_ATA', 'ATA', 'viridis', 'ATA_2baseline_path'),
    ]

    # === Helper: Map context → context_study_tag ===
    def get_context_study_tag(context):
        idx = subject['ASL_CONTEXT'].index(context)
        return subject['context_study_tags'][idx]

    # === Helper: Save NIfTI + DICOM ===
    def save_nifti_and_dicom(context, fields, allow_dicom=True):
        context_study_tag = get_context_study_tag(context)
        
        for field, range_tag, type_tag, _, output_path in fields:
            data = subject[context].get(field) if field != 'CVR_smth' else subject['CVR_smth']
            path = subject[context].get(output_path) if field != 'CVR_smth' else subject['output_CVR_path']
            template = subject[context]['templateNIFTI_path'] if field != 'CVR_smth' else subject['baseline']['templateNIFTI_path']

            if data is not None and path:
                save_data_nifti(data, path, template, 1, None, None)
                
                if allow_dicom:
                    dcm_key = f'templateDCM_{type_tag}_path'
                    dcm_template = subject[context].get(dcm_key, None)
                    dcm_outputdir = subject['DICOMoutputdir']
                    
                    if dcm_template and os.path.exists(dcm_template):
                        save_data_dicom(data,
                                        dcm_template,
                                        dcm_outputdir,
                                        f'ASL {context_study_tag} {type_tag}',  # Here use clinical tag (eg 'ASL preACZ CBF') for SeriesDescription
                                        subject[range_tag],
                                        type_tag
                        )
                    else:
                        logging.warning(f"Skipped DICOM export for type '{type_tag}' ({context_study_tag}): no valid DICOM template found.")
    
    # === Helper: Save PNG ===
    def save_png(context, fields):
        context_study_tag = get_context_study_tag(context)
        for field, range_tag, type_tag, cmap, _ in fields:
            data = subject[context].get(field) if field != 'CVR_smth' else subject['CVR_smth']            
            if data is not None:
                if field == 'CVR_smth':
                    png_name = 'ASL_CVR'
                    title_name = 'CVR'                  
                else:
                    png_name = f'ASL_{context}_{context_study_tag}_{type_tag}'  # e.g., ASL_baseline_preACZ_CBF
                    title_name = f"ASL {context_study_tag} {type_tag}"  # e.g., ASL preACZ CBF
                
                save_figure_to_png(data, subject['nanmask_combined'], subject[range_tag],
                                    subject['RESULTSdir'], png_name, title_name, type_tag, cmap)                
                dcm_template = subject[context][f"templateDCM_{type_tag}_path"]
                save_png_to_dicom(
                    os.path.join(subject['RESULTSdir'], f"{png_name}.png"),
                    os.path.join(subject['RESULTSdir'], f"{png_name}.dcm"), 
                    title_name,
                    dcm_template
                )
                logging.info(f"Saved PNG: {png_name}.png")

    # === Execute saves ===    
    #   Save NIfTI and DICOM for baseline and stimulus
    save_nifti_and_dicom('baseline', fields_main)
    save_nifti_and_dicom('baseline', fields_cvr)
    save_nifti_and_dicom('stimulus', fields_main)
    save_nifti_and_dicom('stimulus', fields_2baseline, allow_dicom=False)

    # Save PNGs for baseline and stimulus
    save_png('baseline', fields_main + fields_cvr)
    save_png('stimulus', fields_2baseline)
        
    # === Save PNGs as DICOM with specific context and type_tag ===    
    # Save PNG as DICOM for baseline CVR
    type_tag = fields_cvr[2] # 'CVR' for baseline
    png_name = f'ASL_{type_tag}'
    input_png_path =  os.path.join(subject['RESULTSdir'], f"{png_name}.png")
    output_dcm_path = os.path.join(subject['RESULTSdir'], f"{png_name}.dcm")
    dcm_template = subject['baseline'][f"templateDCM_{type_tag}_path"]
    series_description = f"ASL {type_tag}" 
    instancenumber = 1  # first instance number for CVR
    save_png_to_dicom(input_png_path, output_dcm_path, series_description, instancenumber, dcm_template)

    # Save PNGs as DICOM for each context and type_tag
    for _, _, type_tag, _, _ in fields_main:
        for context in subject['ASL_CONTEXT']:
            context_study_tag = get_context_study_tag(context)
            png_name = f'ASL_{context}_{context_study_tag}_{type_tag}'
            input_png_path =  os.path.join(subject['RESULTSdir'], f"{png_name}.png")
            output_dcm_path = os.path.join(subject['RESULTSdir'], f"{png_name}.dcm")
            dcm_template = subject[context][f"templateDCM_{type_tag}_path"]
            series_description = f"ASL {context_study_tag} {type_tag}"                  
            instancenumber += 1
            save_png_to_dicom(input_png_path, output_dcm_path, series_description, instancenumber, dcm_template)

    # === Final log ===
    logging.info(f"Results complete: PACS-ready DICOMS, NIFTI, .png's saved for subject {subject['SUBJECTdir']}")
