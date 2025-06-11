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
from clinical_asl_pipeline.asl_t1_from_m0 import asl_t1_from_m0
from clinical_asl_pipeline.asl_qasl_analysis import asl_qasl_analysis
from clinical_asl_pipeline.asl_motioncorrection_ants import asl_motioncorrection_ants
from clinical_asl_pipeline.asl_registration_prepostACZ_ANTS import asl_registration_prepostACZ_ANTS
from clinical_asl_pipeline.asl_save_results_cbfaatcvr import asl_save_results_cbfaatcvr
from clinical_asl_pipeline.utils.littleutils import append_mc
from clinical_asl_pipeline.utils.run_bet_mask import run_bet_mask

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Default processing parameters
DEFAULT_PARAMETERS = {
    'tau': 2,                # Label duration in [s]
    'N_BS': 4,               # Number of background suppression pulses
    'labeleff': 0.85,        # Labeling efficiency (default 0.85)
    'lambda': 0.9,           # Blood-brain partition coefficient (mL/g)
    'T1t': 1.3,              # T1 of gray matter tissue in [s]
    'T1b': 1.65,             # T1 of arterial blood in [s]
    'FWHM': 6,               # Full-width at half-maximum for CVR smoothing [mm]
    'range_cbf': [0, 100],   # Display/analysis range for CBF [ml/100g/min]
    'range_cvr': [-50, 50],  # Display/analysis range for CVR [delta ml/100g/min], cerebrovascular reactivity map upon acetazolamide (DIAMO)
    'range_AAT': [0, 3.5],   # Display/analysis range for AAT [s], arterial arrival time map
    'range_ATA': [0, 125],   # Display/analysis range for ATA [m], arterial transit artefacts map
    'inference_method': 'ssvb', # Inference method for quantification ('ssvb' or 'basil')
}

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
    #
    # Adds keys:
    # - preACZ and postACZ PLD label-control paths
    # - preACZ and postACZ M0 paths
    # - preACZ and postACZ T1fromM0 paths
    # - preACZ and postACZ masks
    # - preACZ and postACZ CBF/AAT/ATA paths
    
    # Initialize phase_tags in subject dictionary
    for phase_tag in ['preACZ', 'postACZ']:
        subject[phase_tag] = {}
    
    # define DICOM phase tags and type tags
    subject['dicom_typetags_by_phasetag'] = {
        'preACZ':  ['CBF', 'AAT', 'ATA', 'CVR'],  # CVR only exists in preACZ space
        'postACZ': ['CBF', 'AAT', 'ATA']  
    }

    # Input data
    for phase_tag in ['preACZ', 'postACZ']:
        # PLD label-control input paths
        subject[phase_tag]['PLDall_labelcontrol_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_allPLD_label1label2.nii.gz')
        subject[phase_tag]['PLD2tolast_labelcontrol_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_2tolastPLD_label1label2.nii.gz')
        subject[phase_tag]['PLD1to2_labelcontrol_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_1to2PLD_label1label2.nii.gz')

        # Mask, M0 and T1fromM0 paths
        subject[phase_tag]['mask_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_M0_brain_mask.nii.gz')
        subject[phase_tag]['M0_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_M0.nii.gz')
        subject[phase_tag]['T1fromM0_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_T1fromM0.nii.gz')

        # QASL ASL derived quantification output paths after QASL analysis for CBF/AAT/ATA
        subject[phase_tag]['QASL_CBF_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_QASL_2tolastPLD_forCBF/output/native/calib_voxelwise/perfusion.nii.gz')
        subject[phase_tag]['QASL_AAT_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_QASL_allPLD_forAAT/output/native/arrival.nii.gz')
        subject[phase_tag]['QASL_ATA_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_QASL_1to2PLD_forATA/output/native/calib_voxelwise/perfusion.nii.gz')

        # ASL derived quantification output paths in ASLdir for CBF/AAT/ATA/CVR
        subject[phase_tag]['output_CBF_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_CBF.nii.gz')
        subject[phase_tag]['output_AAT_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_AAT.nii.gz')
        subject[phase_tag]['output_ATA_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_ATA.nii.gz')
        subject['output_CVR_path'] = os.path.join(subject['ASLdir'], 'CVR.nii.gz') #  in parent subject dictionary as it is considers both pre- and postACZ data

        # Registered postACZ CBF/AAT/ATA/M0/T1fromM0/mask → preACZ space — only for postACZ
        if phase_tag == 'postACZ':
            subject[phase_tag]['CBF_2preACZ_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_CBF_2preACZ.nii.gz')
            subject[phase_tag]['AAT_2preACZ_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_AAT_2preACZ.nii.gz')
            subject[phase_tag]['ATA_2preACZ_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_ATA_2preACZ.nii.gz')
            subject[phase_tag]['M0_2preACZ_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_M0_2preACZ.nii.gz')
            subject[phase_tag]['T1fromM0_2preACZ_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_T1fromM0_2preACZ.nii.gz')
            subject[phase_tag]['mask_2preACZ_path'] = os.path.join(subject['ASLdir'], f'{phase_tag}_M0_brain_mask_2preACZ.nii.gz')

    logging.info("Input and derived ASL file paths prepared.")
    return subject

def get_latest_source_data(dicomdir, niftidir, phase_tag):
    # Find the latest SOURCE_ASL DICOM and NIfTI files for a given tag (e.g., 'preACZ', 'postACZ').
    # Parameters:
    #     dicomdir (str): Path to the DICOM directory 
    #     niftidir (str): Path to the NIfTI directory
    #     phase_tag (str): Tag to identify phase ('preACZ' or 'postACZ')
    #
    # Returns:
    #     tuple: (nifti_full_path, dicom_full_path)

    # Find DICOM file names (ending in '2', e.g., series file names)
    dicom_files = sorted([f for f in os.listdir(dicomdir) if 'SOURCE_ASL' in f and phase_tag in f and f.endswith('2')])
    if len(dicom_files) > 1:
        warnings.warn(f'Multiple SOURCE_ASL DICOM entries found for "{phase_tag}". Using latest.')
    if not dicom_files:
        raise FileNotFoundError(f'No SOURCE_ASL DICOM entry found for phase_tag "{phase_tag}" in {dicomdir}')

    dicom_path = os.path.join(dicomdir, dicom_files[-1])

    # Find NIFTI file names (ending in '2.nii.gz', e.g., series file names)
    files_nifti = sorted([f for f in os.listdir(niftidir) if 'SOURCE_ASL' in f and phase_tag in f and f.endswith('2.nii.gz')])
    if len(files_nifti) > 1:
        warnings.warn(f'Multiple SOURCE_ASL NIFTI files found for "{phase_tag}". Using latest.')
    if not files_nifti:
        raise FileNotFoundError(f'No SOURCE_ASL NIFTI files found for phase_tag "{phase_tag}" in {niftidir}')

    nifti_path = os.path.join(niftidir, files_nifti[-1])

    return nifti_path, dicom_path

def find_template_dicom(files, type_tag, phase_tag): # Helper function to find a template DICOM file based on type and phase tags
    return next(f for f in files if 'sWIP' in f and type_tag in f and phase_tag in f)


def mri_diamox_umcu_clinicalasl_cvr_imager(inputdir, outputdir):
    # Main function to run the Clinical ASL pipeline for a subject.
    # Parameters:
    #     inputdir (str): Path to the input directory containing extracted PACS DICOM files.
    #     outputdir (str): Path to the output directory where ASL derived images and generated DICOMS will be saved.
    # Initialize subject dictionary with input and output directories

    logging.info("Starting Clinical ASL pipeline for subject...")

    subject = {}
    subject['DICOMinputdir'] = inputdir
    subject['SUBJECTdir'] = outputdir 
    
    ##### Step 0: Prepare subject dictionary with default parameters and paths
    # Add default parameters to subject dictionary, such as scan parameters, FWHM, ranges, etc.
    subject.update(DEFAULT_PARAMETERS)

    # Prepare output folders
    subject = prepare_subject_paths(subject)

    # Prepare file paths
    subject = prepare_input_output_paths(subject)

    ###### Step 1: Convert DICOM to NIFTI, move input PACS DICOMSinputdir to DICOMsubjectdir for further processing
    subject = asl_convert_dicom_to_nifti(subject)

    ###### Step 2: Get SOURCE and DICOM NIFTI files using phase tags 
    for phase_tag in ['preACZ', 'postACZ']:
        nifti_file, dicom_file = get_latest_source_data(subject['DICOMsubjectdir'], subject['NIFTIdir'], phase_tag)
        nifti_full_path = os.path.join(subject['NIFTIdir'], nifti_file)
        dicom_full_path = os.path.join(subject['DICOMsubjectdir'], dicom_file)
        subject[phase_tag]['sourceNIFTI_path'] = nifti_full_path
        subject[phase_tag]['templateNII_path'] = nifti_full_path
        subject[phase_tag]['sourceDCM_path']   = dicom_full_path

    ###### Step 3: Locate template DICOMs - define DICOM phase tags and type tags
    for phase_tag, type_tags in subject['dicom_typetags_by_phasetag'].items():
        for type_tag in type_tags:
            subject[phase_tag][f'templateDCM_{type_tag}_path'] = find_template_dicom(os.listdir(subject['DICOMsubjectdir']), type_tag, phase_tag)

    ###### Step 4: DICOM scanparameter extraction 
    subject = asl_extract_params_dicom(subject, subject['preACZ']['sourceDCM_path'], phase_tag='preACZ')
    subject = asl_extract_params_dicom(subject, subject['postACZ']['sourceDCM_path'], phase_tag='postACZ')

    ###### Step 5: Compute Look Locker correction on Mxy per PLD   
    subject = asl_look_locker_correction(subject, phase_tag='preACZ')
    subject = asl_look_locker_correction(subject, phase_tag='postACZ')

    ###### Step 6: Split/control-label, save NIFTI
    subject = asl_prepare_asl_data(subject, subject['preACZ']['sourceNIFTI_path'], phase_tag='preACZ')
    subject = asl_prepare_asl_data(subject, subject['postACZ']['sourceNIFTI_path'], phase_tag='postACZ')

    ###### Step 7: Brain extraction on M0 images using HD-BET CLI
    subject['preACZ']['mask'], subject['preACZ']['nanmask']  = run_bet_mask(subject['preACZ']['M0_path'], subject['preACZ']['mask_path'])
    subject['postACZ']['mask'], subject['postACZ']['nanmask'] = run_bet_mask(subject['postACZ']['M0_path'], subject['postACZ']['mask_path'])
    
    ###### Step 8:# compute T1 from M0
    subject = asl_t1_from_m0(subject, phase_tag='preACZ') 
    subject = asl_t1_from_m0(subject, phase_tag='postACZ')

    #for phase_tag in ['preACZ', 'postACZ']:    
    ###### Step 9: Motion Correction on label and control images
        #asl_motioncorrection_ants(subject[phase_tag]['PLDall_labelcontrol_path'], subject[phase_tag]['M0_path'], append_mc(subject[phase_tag]['PLDall_labelcontrol_path']))
        #asl_motioncorrection_ants(subject[phase_tag]['PLD2tolast_labelcontrol_path'], subject[phase_tag]['M0_path'], append_mc(subject[phase_tag]['PLD2tolast_labelcontrol_path']))
        #asl_motioncorrection_ants(subject[phase_tag]['PLD1to2_labelcontrol_path'], subject[phase_tag]['M0_path'], append_mc(subject[phase_tag]['PLD1to2_labelcontrol_path']))        

    ###### Step 10: ASL Quantification analysis
        # all PLD for AAT (arterial arrival time map)
        #asl_qasl_analysis(subject, append_mc(subject[phase_tag]['PLDall_labelcontrol_path']), subject[phase_tag]['M0_path'], subject[phase_tag]['mask_path'] , os.path.join(subject['ASLdir'], f'{phase_tag}_QASL_allPLD_forAAT'), subject['PLDS'][0:], subject['inference_method'])
       
        # 2-to-last PLD for CBF map
        #asl_qasl_analysis(subject, append_mc(subject[phase_tag]['PLD2tolast_labelcontrol_path']), subject[phase_tag]['M0_path'], subject[phase_tag]['mask_path'] , os.path.join(subject['ASLdir'], f'{phase_tag}_QASL_2tolastPLD_forCBF'), subject['PLDS'][1:], subject['inference_method'])
       
        # 1to2 PLDs for ATA map ->  then do no fit for the arterial component 'artoff'
        #asl_qasl_analysis(subject, append_mcsubject[phase_tag]['PLD1to2_labelcontrol_path']), subject[phase_tag]['M0_path'], subject[phase_tag]['mask_path'] , os.path.join(subject['ASLdir'], f'{phase_tag}_QASL_1to2PLD_forATA'), subject['PLDS'][0:2], subject['inference_method'], 'artoff')
    
    ###### Step 11: register post-ACZ ASL data to pre-ACZ ASL data using Elastix 
    asl_registration_prepostACZ_ANTS(subject)

    ###### Step 12: Generate CBF/AAT/ATA/CVR results (nifti, dicom PACS, .pngs) for pre- and postACZ, including registration of postACZ to preACZ as reference data, and target for computed CVR map
    subject = asl_save_results_cbfaatcvr(subject)

    logging.info("ASL processing pipeline completed successfully.")
    # Return the subject dictionary with all paths and results
    return subject

  