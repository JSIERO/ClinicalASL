"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Main pipeline script for processing ASL MRI data.
Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    This script implements the main pipeline for processing ASL data, including DICOM to NIfTI conversion,
    parameter extraction, Look-Locker correction, data preparation, brain extraction, motion correction,
    quantification, registration, and saving of CBF, AAT, ATA, and CVR results for a subject.
    Results are saved as NIfTI, PNG, and DICOM files for further analysis and PACS export.

License: BSD 3-Clause License
"""

import os
import logging
import warnings
from clinical_asl_pipeline.asl_convert_dicom_to_nifti import asl_convert_dicom_to_nifti
from clinical_asl_pipeline.asl_extract_params_dicom import asl_extract_params_dicom
from clinical_asl_pipeline.asl_look_locker_correction import asl_look_locker_correction
from clinical_asl_pipeline.asl_prepare_asl_data import asl_prepare_asl_data
from clinical_asl_pipeline.asl_motion_correction import asl_motion_correction
from clinical_asl_pipeline.asl_outlier_removal import asl_outlier_removal
from clinical_asl_pipeline.asl_t1_from_m0 import asl_t1_from_m0
from clinical_asl_pipeline.asl_qasl_analysis import asl_qasl_analysis
from clinical_asl_pipeline.asl_registration_stimulus_to_baseline import asl_registration_stimulus_to_baseline
from clinical_asl_pipeline.asl_save_results_cbfaatcvr import asl_save_results_cbfaatcvr
from clinical_asl_pipeline.utils.run_bet_mask import run_bet_mask

def prepare_subject_paths(subject):
    # Prepare output folder structure for subject.
    #
    # Creates the following subfolders:
    # - NIFTI
    # - ASL
    # - FIGURE_RESULTS
    # - (reuses SUBJECTdir as DICOMoutputdir)
    subject['DICOMoutputdir'] = subject['SUBJECTdir'] #     
    subject['DICOMsubjectdir'] = os.path.join(subject['SUBJECTdir'], 'DICOMORIG') # folder for original DICOMS in subject output fodler
    subject['NIFTIdir'] = os.path.join(subject['SUBJECTdir'], 'NIFTI')
    subject['ASLdir'] = os.path.join(subject['SUBJECTdir'], 'ASL')
    subject['RESULTSdir'] = os.path.join(subject['SUBJECTdir'], 'FIGURE_RESULTS')

    for path in [
        subject['DICOMoutputdir'],
        subject['DICOMsubjectdir'],
        subject['NIFTIdir'],
        subject['ASLdir'],
        subject['RESULTSdir']
    ]:
        os.makedirs(path, exist_ok=True)

    logging.info("Output folder structure created.")
    return subject

def prepare_input_output_paths(subject):
    # Prepare standard input and derived ASL file paths in subject dictionary.
    # Initialize context in subject dictionary
    for context in subject['ASL_CONTEXT']:
        subject[context] = {}

    # define DICOM context tags and type tags
    subject['dicom_typetags_by_context'] = {
        'baseline':  ['CBF', 'AAT', 'ATA', 'CVR'],
        'stimulus': ['CBF', 'AAT', 'ATA']
    }

    # Input data
    for context in subject['ASL_CONTEXT']:
        subject[context]['PLDall_controllabel_path'] = os.path.join(subject['ASLdir'], f'{context}_allPLD_controllabel.nii.gz')
        subject[context]['PLD2tolast_controllabel_path'] = os.path.join(subject['ASLdir'], f'{context}_2tolastPLD_controllabel.nii.gz')
        subject[context]['PLD1to2_controllabel_path'] = os.path.join(subject['ASLdir'], f'{context}_1to2PLD_controllabel.nii.gz')

        subject[context]['mask_path'] = os.path.join(subject['ASLdir'], f'{context}_M0_brain_mask.nii.gz')
        subject[context]['M0_path'] = os.path.join(subject['ASLdir'], f'{context}_M0.nii.gz')
        subject[context]['T1fromM0_path'] = os.path.join(subject['ASLdir'], f'{context}_T1fromM0.nii.gz')

        subject[context]['QASL_CBF_path'] = os.path.join(subject['ASLdir'], f'{context}_QASL_2tolastPLD_forCBF/output/native/calib_voxelwise/perfusion.nii.gz')
        subject[context]['QASL_AAT_path'] = os.path.join(subject['ASLdir'], f'{context}_QASL_allPLD_forAAT/output/native/arrival.nii.gz')
        subject[context]['QASL_ATA_path'] = os.path.join(subject['ASLdir'], f'{context}_QASL_1to2PLD_forATA/output/native/calib_voxelwise/perfusion.nii.gz')

        subject[context]['output_CBF_path'] = os.path.join(subject['ASLdir'], f'{context}_CBF.nii.gz')
        subject[context]['output_AAT_path'] = os.path.join(subject['ASLdir'], f'{context}_AAT.nii.gz')
        subject[context]['output_ATA_path'] = os.path.join(subject['ASLdir'], f'{context}_ATA.nii.gz')

        subject['output_CVR_path'] = os.path.join(subject['ASLdir'], 'CVR.nii.gz')

        if context == 'stimulus':
            subject[context]['CBF_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_CBF_2baseline.nii.gz')
            subject[context]['AAT_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_AAT_2baseline.nii.gz')
            subject[context]['ATA_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_ATA_2baseline.nii.gz')
            subject[context]['M0_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_M0_2baseline.nii.gz')
            subject[context]['T1fromM0_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_T1fromM0_2baseline.nii.gz')
            subject[context]['mask_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_M0_brain_mask_2baseline.nii.gz')

    logging.info("Input and derived ASL file paths prepared.")
    return subject

def get_latest_source_data(dicomdir, niftidir, context_study_tag):
    # Find the latest SOURCE_ASL DICOM and NIfTI files for a given tag (e.g., 'baseline' or 'stimulus')
    # Parameters:
    #     dicomdir (str): Path to the DICOM directory 
    #     niftidir (str): Path to the NIfTI directory
    #     context_study_tag (str): Study-specific context tag to identify filenames for context 'baseline' or 'stimulus', defined at ANALYSIS_PARAMETERS
    # Returns:
    #     tuple: (nifti_full_path, dicom_full_path)

    # Find DICOM file names (ending in '2', e.g., series file names)
    dicom_files = sorted([f for f in os.listdir(dicomdir) if 'SOURCE_ASL' in f and context_study_tag in f and f.endswith('2')])
    if len(dicom_files) > 1:
        warnings.warn(f'Multiple SOURCE ASL DICOM entries found for "{context_study_tag}". Using latest.')
    if not dicom_files:
        raise FileNotFoundError(f'No SOURCE ASL DICOM entry found for context "{context_study_tag}" in {dicomdir}')

    dicom_path = os.path.join(dicomdir, dicom_files[-1])

    # Find NIFTI file names (ending in '2.nii.gz', e.g., series file names)
    files_nifti = sorted([f for f in os.listdir(niftidir) if 'SOURCE_ASL' in f and context_study_tag in f and f.endswith('2.nii.gz')])
    if len(files_nifti) > 1:
        warnings.warn(f'Multiple SOURCE ASL NIFTI files found for "{context_study_tag}". Using latest.')
    if not files_nifti:
        raise FileNotFoundError(f'No SOURCE ASL NIFTI files found for context "{context_study_tag}" in {niftidir}')

    nifti_path = os.path.join(niftidir, files_nifti[-1])

    return nifti_path, dicom_path

def find_template_dicom(files, type_tag, context_study_tag): 
    # Helper function to find a template DICOM file based on type and context tags
    return next(f for f in files if 'sWIP' in f and type_tag in f and context_study_tag in f)

def mri_diamox_umcu_clinicalasl_cvr(inputdir, outputdir, ANALYSIS_PARAMETERS):
    # Main function to run the Clinical ASL pipeline for a subject.
    # Parameters:
    #     inputdir (str): Path to the input directory containing extracted PACS DICOM files.
    #     outputdir (str): Path to the output directory where ASL derived images and generated DICOMS will be saved.
    # Initialize subject dictionary with input and output directories

    subject = {}
    subject['DICOMinputdir'] = inputdir
    subject['SUBJECTdir'] = outputdir 

    ##### Step 1: Prepare subject dictionary with default parameters and paths
    # Add default parameters to subject dictionary, such as context tags for baseline/stimulus scans, scan parameters, FWHM, ranges, etc.
    subject.update(ANALYSIS_PARAMETERS)

    # Prepare output folders
    subject = prepare_subject_paths(subject)

    # Prepare file paths
    subject = prepare_input_output_paths(subject)

    ###### Step 2: Convert DICOM to NIFTI, move input PACS DICOMSinputdir to DICOMsubjectdir for further processing
    subject = asl_convert_dicom_to_nifti(subject)

    ###### Step 3-11: Unified loop for each ASL context tag: 'baseline', 'stimulus'
    for i, context in enumerate(subject['ASL_CONTEXT']):
        context_study_tag = subject['context_study_tags'][i]
        logging.info(f"=== Processing context '{context}' (tag: '{context_study_tag}') ===")

    ###### Step 3: Get SOURCE and DICOM NIFTI files
        nifti_file, dicom_file = get_latest_source_data(subject['DICOMsubjectdir'], subject['NIFTIdir'], context_study_tag=context_study_tag)
        subject[context]['sourceNIFTI_path'] = os.path.join(subject['NIFTIdir'], nifti_file)
        subject[context]['templateNII_path'] = os.path.join(subject['NIFTIdir'], nifti_file)
        subject[context]['sourceDCM_path']   = os.path.join(subject['DICOMsubjectdir'], dicom_file)

    ###### Step 4: Locate template DICOMs
        context_study_tag = subject['context_study_tags'][i]
        for type_tag in subject['dicom_typetags_by_context'][context]:
            subject[context][f'templateDCM_{type_tag}_path'] = find_template_dicom(os.listdir(subject['DICOMsubjectdir']), type_tag, context_study_tag)

    ###### Step 5: DICOM scanparameter extraction
        subject = asl_extract_params_dicom(subject, subject[context]['sourceDCM_path'], context_tag=context)

    ###### Step 6: Look Locker correction
        subject = asl_look_locker_correction(subject, context_tag=context)

    ###### Step 7: Interleave control-label, save to NIFTI
        subject = asl_prepare_asl_data(subject, subject[context]['sourceNIFTI_path'], context_tag=context)

    ###### Step 8: Brain extraction on M0 using HD-BET CLI
        subject[context]['mask'], subject[context]['nanmask'] = run_bet_mask(subject[context]['M0_path'], subject[context]['mask_path'], device = subject['device'], extradata_path = subject[context]['PLDall_controllabel_path'])
    
        subject = asl_motion_correction(subject, context_tag=context)

    ###### Step 8: Outlier timepoint rejection: 2.5 x std + mean CBF (deltaM) 
        subject = asl_outlier_removal(subject, context_tag=context, usermask=None)

    ###### Step 9: Compute T1 from M0
        subject = asl_t1_from_m0(subject, context_tag=context)

    ###### Step 10: ASL Quantification analysis
        context_data = subject[context]
        # # all PLD for AAT (arterial arrival time map)
        asl_qasl_analysis(context_data, ANALYSIS_PARAMETERS, 
                        context_data['PLDall_controllabel_path'], 
                        context_data['M0_path'], 
                        context_data['mask_path'], 
                        os.path.join(subject['ASLdir'], f'{context}_QASL_allPLD_forAAT'),       # output folder name QASL
                        context_data['PLDS'][0:], 
                        subject['inference_method']
                        )
        # 2-to-last PLD for CBF map
        asl_qasl_analysis(context_data, ANALYSIS_PARAMETERS, 
                        context_data['PLD2tolast_controllabel_path'], 
                        context_data['M0_path'], 
                        context_data['mask_path'], 
                        os.path.join(subject['ASLdir'], f'{context}_QASL_2tolastPLD_forCBF'),   # output folder name QASL
                        subject[context]['PLDS'][1:], 
                        subject['inference_method']
                        )
        # 1to2 PLDs for ATA map ->  then do no fit for the arterial component 'artoff'
        asl_qasl_analysis(context_data, ANALYSIS_PARAMETERS, 
                        context_data['PLD1to2_controllabel_path'], 
                        context_data['M0_path'], 
                        context_data['mask_path'], 
                        os.path.join(subject['ASLdir'], f'{context}_QASL_1to2PLD_forATA'),      # output folder name QASL
                        context_data['PLDS'][0:2], 
                        subject['inference_method'],
                        'artoff'
                        )

    ###### Step 11: register post-ACZ ASL data to pre-ACZ ASL data using Elastix 
    asl_registration_stimulus_to_baseline(subject)

    ###### Step 12: Generate CBF/AAT/ATA/CVR results (nifti, dicom PACS, .pngs) for pre- and postACZ, including registration of postACZ to preACZ as reference data, and target for computed CVR map
    asl_save_results_cbfaatcvr(subject)

    logging.info("ASL processing pipeline completed successfully.")
