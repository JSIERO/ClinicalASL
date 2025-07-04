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
import pydicom
import fnmatch
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
    # - (uses SUBJECTdir as DICOMoutputdir)
    subject['DICOMoutputdir'] = subject['SUBJECTdir'] # folder were the derived CBF, AAT, CVR, ATA DICOMS will be saved
    subject['DICOMsubjectdir'] = os.path.join(subject['SUBJECTdir'], 'DICOMORIG') # folder location of DICOMS for processing (copied and renamed from PACS or scanner)
    subject['DICOMorigdir'] = os.path.join(subject['DICOMsubjectdir'], 'ORIG') # folder location of exact copy of original DICOMS as copied from pACS or Philips scanner
    subject['NIFTIdir'] = os.path.join(subject['SUBJECTdir'], 'NIFTI')
    subject['ASLdir'] = os.path.join(subject['SUBJECTdir'], 'ASL')
    subject['RESULTSdir'] = os.path.join(subject['SUBJECTdir'], 'FIGURE_RESULTS')

    for path in [
        subject['DICOMoutputdir'],
        subject['DICOMsubjectdir'],
        subject['DICOMorigdir'] ,
        subject['NIFTIdir'],
        subject['ASLdir'],
        subject['RESULTSdir']
    ]:
        os.makedirs(path, exist_ok=True)

    logging.info("Output folder structure created.")
    return subject

def prepare_input_output_paths(subject):
    # Prepare standard input and derived ASL file paths in subject dictionary.
    # Initialize context in subject dictionary; default 'baseline', 'stimulus'
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

        subject[context]['mask_path'] = os.path.join(subject['ASLdir'], f'{context}_M0_brain_mask.nii.gz')
        subject[context]['M0_path'] = os.path.join(subject['ASLdir'], f'{context}_M0.nii.gz')

        subject[context]['QASL_CBF_path'] = os.path.join(subject['ASLdir'], f'{context}_QASL/output/native/calib_voxelwise/perfusion.nii.gz')
        subject[context]['QASL_AAT_path'] = os.path.join(subject['ASLdir'], f'{context}_QASL/output/native/arrival.nii.gz')

        subject[context]['output_CBF_path'] = os.path.join(subject['ASLdir'], f'{context}_CBF.nii.gz')
        subject[context]['output_AAT_path'] = os.path.join(subject['ASLdir'], f'{context}_AAT.nii.gz')

        subject['output_CVR_path'] = os.path.join(subject['ASLdir'], 'CVR.nii.gz')

        if context == 'stimulus':
            subject[context]['CBF_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_CBF_2baseline.nii.gz')
            subject[context]['AAT_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_AAT_2baseline.nii.gz')
            subject[context]['M0_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_M0_2baseline.nii.gz')
            subject[context]['mask_2baseline_path'] = os.path.join(subject['ASLdir'], f'{context}_M0_brain_mask_2baseline.nii.gz')

    logging.info("Input and derived ASL file paths prepared.")
    return subject

def get_latest_source_data(subject, context_study_tag, context_tag):
    # Find the latest matching ASL DICOM and NIfTI files for a given tag (e.g., 'baseline' or 'stimulus')
    # based on SeriesDescription patterns.
    #
    # This function:
    # - Searches the DICOM and NIfTI directories for files containing a user-defined context tag 
    #   (e.g., 'baseline', 'stimulus') and matching a specified SeriesDescription pattern.
    # - By default, the SeriesDescription patterns matched are ['*SOURCE*ASL*'] (case-insensitive).
    # - Picks the latest file (based on sorted filename) when multiple matches are found.
    #
    # Parameters:
    #     dicomdir (str): Path to the DICOM directory to search.
    #     niftidir (str): Path to the NIfTI directory to search.
    #     context_study_tag (str): Study-specific context tag for selecting appropriate files.
    #     context_tag: string indicating the ASL context tag for the subject's data, e.g., 'baseline', 'stimulus', etc.
    #     subject (dict): Subject dictionary that may contain optional
    #         'dicomseries_description_patterns' key with custom matching patterns.
    #
    # Returns:
    #     tuple: (nifti_full_path, dicom_full_path) corresponding to the latest matching entries.

    dicomdir = subject['DICOMsubjectdir']
    niftidir = subject['NIFTIdir']
    context_data = subject[context_tag]

    series_patterns = subject.get('dicomseries_description_patterns', ['*SOURCE*vTR*', '*SOURCE*M0*'])

    # Filter DICOMs
    if subject['is_singleframe']:
        series_files = {}
        for fname in os.listdir(dicomdir):
            if context_study_tag not in fname or not fnmatch.fnmatch(fname.upper(), series_patterns[0].upper()):
                continue
            fpath = os.path.join(dicomdir, fname)
            try:
                ds = pydicom.dcmread(fpath, stop_before_pixels=True, force=True)
                sn = int(getattr(ds, 'SeriesNumber', -1))
                if sn not in series_files:
                    series_files[sn] = []
                series_files[sn].append(fname)
            except Exception as e:
                continue  # skip unreadable or invalid DICOMs

        if not series_files:
            raise FileNotFoundError(f"No matching DICOM files found for context '{context_study_tag}' in {dicomdir}")

        highest_sn = max(series_files)
        selected_dicom_file = sorted(series_files[highest_sn])[0]  # get the first slice in sorted order

        dicom_path = os.path.join(dicomdir, selected_dicom_file)

    elif subject['is_multiframe']:
        #filter DICOMs by SeriesDescription patterns
        def select_dicom_file(pattern_idx, context_tag, description):
            # Helper to select DICOM file by SeriesNumber for given pattern index and context, for SOURCE ASL AND M0 respectively
            matched_files = sorted([
                fname for fname in os.listdir(dicomdir)
                if fnmatch.fnmatch(fname.upper(), series_patterns[pattern_idx].upper()) and fname.endswith('2')
            ])
            if len(matched_files) > 1:
                logging.warning(f"Multiple {description} DICOM entries found. Sorting according to SeriesNumber, using 1st entry for context_tag=baseline, 2nd entry for context_tag=stimulus.")
                matched_files = sorted(
                    matched_files,
                    key=lambda x: int(
                        pydicom.dcmread(os.path.join(dicomdir, x), stop_before_pixels=True, force=True).SeriesNumber
                    )
                )
                if context_tag == 'baseline':
                    return os.path.join(dicomdir, matched_files[0])
                elif context_tag == 'stimulus':
                    return os.path.join(dicomdir, matched_files[-1])
            elif len(matched_files) == 1:
                logging.warning(f"Only one {description} DICOM entry found for \"{context_study_tag}\". Using it.")
                return os.path.join(dicomdir, matched_files[0])
            else:
                raise FileNotFoundError(f"No {description} DICOM entry found for context \"{context_study_tag}\" in {dicomdir}")

        dicom_path = select_dicom_file(0, context_tag, "SOURCE ASL")
        dicom_m0_path = select_dicom_file(1, context_tag, "SOURCE M0")

        # Filter NIFTIs
        def select_nifti_file(pattern_idx, context_tag, description):
            # Helper to select NIFTI file by SeriesNumber for given pattern index and context, for SOURCE ASL AND M0 respectively
            matched_files = sorted([
            fname for fname in os.listdir(niftidir)
            if fnmatch.fnmatch(fname.upper(), series_patterns[pattern_idx].upper()) and fname.endswith('2.nii.gz')
            ])
            if len(matched_files) > 1:
                logging.warning(f"Multiple {description} NIFTI entries found. Sorting according to SeriesNumber, using 1st entry for context_tag=baseline, 2nd entry for context_tag=stimulus.")
            def extract_series_number(filename):
                import re
                # Extract digits after last underscore before .nii.gz
                match = re.search(r'_([0-9]+)\.nii\.gz$', filename)
                return int(match.group(1)) if match else -1

            matched_files = sorted(matched_files, key=extract_series_number)
            if context_tag == 'baseline':
                return os.path.join(niftidir, matched_files[0])
            elif context_tag == 'stimulus':
                return os.path.join(niftidir, matched_files[-1])
            elif len(matched_files) == 1:
                logging.warning(f"Only one {description} NFITI entry found for \"{context_study_tag}\". Using it.")
                return os.path.join(niftidir, matched_files[0])
            else:
                raise FileNotFoundError(f"No {description} NIFTI entry found for context \"{context_study_tag}\" in {niftidir}")
            
        nifti_path = select_nifti_file(0, context_tag, "SOURCE ASL")
        nifti_m0_path = select_nifti_file(1, context_tag, "SOURCE M0")
        
    context_data['sourceNIFTI_path'] = nifti_path
    context_data['sourceNIFTI_path'] = nifti_path
    context_data['sourceDCM_path']   = dicom_path
    context_data['sourceM0NIFTI_path'] = nifti_m0_path
    context_data['sourceM0DCM_path']   = dicom_m0_path

    logging.info(f"source DICOM file selected: {dicom_path}")
    logging.info(f"source NIFTI file selected: {nifti_path}")
    logging.info(f"source M0 DICOM file selected: {dicom_m0_path}") 
    logging.info(f"source M0 NIFTI file selected: {nifti_m0_path}")

    return subject

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
        context_study_tag = subject['context_study_tags'][i] # e.g 'preACZ' and 'postACZ'

        logging.info(f"===================================================================")
        logging.info(f"=== Processing context '{context}' (tag: '{context_study_tag}') ===")
        logging.info(f"===================================================================")

    ###### Step 3: Get SOURCE and DICOM NIFTI files
        subject = get_latest_source_data(subject, context_study_tag, context_tag=context)

    ###### Step 5: DICOM scanparameter extraction
        subject = asl_extract_params_dicom(subject, context_tag=context)

    ###### Step 7: Interleave control-label, save to NIFTI
        subject = asl_prepare_asl_data(subject, context_tag=context)

    ###### Step 8: Brain extraction on M0 using HD-BET CLI
        subject = run_bet_mask(subject, context_tag=context)
    
        subject = asl_motion_correction(subject, context_tag=context)

    ###### Step 8: Outlier timepoint rejection: 2.5 x std + mean CBF (deltaM) 
        subject = asl_outlier_removal(subject, context_tag=context, usermask=None)

    ###### Step 10: ASL Quantification analysis
        context_data = subject[context]
        asl_qasl_analysis(context_data, ANALYSIS_PARAMETERS, 
                        context_data['PLDall_controllabel_path'], 
                        context_data['M0_path'], 
                        context_data['mask_path'], 
                        os.path.join(subject['ASLdir'], f'{context}_QASL_allPLD_forAAT'),       # output folder name QASL
                        context_data['PLDS'][0:], 
                        subject['inference_method']
                        )


    ###### Step 11: register post-ACZ ASL data to pre-ACZ ASL data using Elastix 
    asl_registration_stimulus_to_baseline(subject)

    ###### Step 12: Generate CBF/AAT/ATA/CVR results (nifti, dicom PACS, .pngs) for pre- and postACZ, including registration of postACZ to preACZ as reference data, and target for computed CVR map
    asl_save_results_cbfaatcvr(subject)

    logging.info("ASL processing pipeline completed successfully.")
