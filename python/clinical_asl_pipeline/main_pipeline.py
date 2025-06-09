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
from clinical_asl_pipeline.asl_convert_dicom_to_nifti import asl_convert_dicom_to_nifti
from clinical_asl_pipeline.asl_extract_params_dicom import asl_extract_params_dicom
from clinical_asl_pipeline.asl_look_locker_correction import asl_look_locker_correction
from clinical_asl_pipeline.asl_prepare_asl_data import asl_prepare_asl_data
from clinical_asl_pipeline.asl_bet_t1_from_m0 import asl_bet_t1_from_m0
from clinical_asl_pipeline.asl_qasl_analysis import asl_qasl_analysis
from clinical_asl_pipeline.asl_motioncorrection_ants import asl_motioncorrection_ants
from clinical_asl_pipeline.asl_registration_prepostACZ import asl_registration_prepostACZ
from clinical_asl_pipeline.asl_save_results_cbfaatcvr import asl_save_results_cbfaatcvr

def mri_diamox_umcu_clinicalasl_cvr_imager(inputdir, outputdir):
    subject = {}
    subject['DICOMdir'] = inputdir # input folder for extracted PACS DICOM
    subject['SUBJECTdir'] = outputdir # inpput folder for results and generated DICOMS of ALS derived imaged for PACS
    
    # Set registration (Elastix) settings
    subject['elastix_parameter_file'] = 'Par0001rigid_6DOF_MI_NIFTIGZ.txt'

    # Default parameters
    subject.update({
        'tau': 2,               # in [s]
        'N_BS': 4,              # in [number of BS]
        'labeleff': 0.85,       # in [0.85=85%]
        'lambda': 0.9,          # in [g/ml]
        'T1t': 1.3,             # in [s]
        'T1b': 1.65,            # in [s]
        'FWHM': 6,              # in [mm]
        'FWHM_M0': 5,           # in [mm]
        'range_cbf': [0, 100],  # in [ml/100g/min]
        'range_cvr': [-50, 50], # in [ml/100g/min]
        'range_AAT': [0, 3.5],  # in [s]
        'range_ATA': [0, 125],  # in [ml/100g/min]
        'inference_method': 'ssvb', # method for ASL quantification, choices: ssvb, basil
    })

    # Set Folder paths
    subject['DICOMRESULTSdir'] = subject['SUBJECTdir'] # root outputfolder for generated DICOMS derived images for PACS
    subject['NIFTIdir'] = os.path.join(subject['SUBJECTdir'], 'NIFTI')
    subject['ASLdir'] = os.path.join(subject['SUBJECTdir'], 'ASL')
    subject['RESULTSdir'] = os.path.join(subject['SUBJECTdir'], 'FIGURE_RESULTS')

    for path in [subject['DICOMRESULTSdir'], subject['NIFTIdir'], subject['ASLdir'], subject['RESULTSdir']]:
        os.makedirs(path, exist_ok=True)
        
    # Set necessary file paths
    # Input data
    subject['preACZ_PLDall_labelcontrol_path'] = os.path.join(subject['ASLdir'], 'preACZ_allPLD_label1label2.nii.gz') 
    subject['preACZ_PLD2tolast_labelcontrol_path'] = os.path.join(subject['ASLdir'], 'preACZ_2tolastPLD_label1label2.nii.gz')
    subject['preACZ_PLD1to2_labelcontrol_path'] = os.path.join(subject['ASLdir'], 'preACZ_1to2PLD_label1label2.nii.gz')

    subject['postACZ_PLDall_labelcontrol_path'] = os.path.join(subject['ASLdir'], 'postACZ_allPLD_label1label2.nii.gz')
    subject['postACZ_PLD2tolast_labelcontrol_path'] = os.path.join(subject['ASLdir'], 'postACZ_2tolastPLD_label1label2.nii.gz')
    subject['postACZ_PLD1to2_labelcontrol_path'] = os.path.join(subject['ASLdir'], 'postACZ_1to2PLD_label1label2.nii.gz')
    
    subject['preACZ_M0_path'] = os.path.join(subject['ASLdir'], 'preACZ_M0.nii.gz')
    subject['postACZ_M0_path'] = os.path.join(subject['ASLdir'], 'postACZ_M0.nii.gz')
    subject['postACZ_M0_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_M0_2preACZ.nii.gz')    

    # Derived data
    subject['preACZ_T1fromM0_path'] = os.path.join(subject['ASLdir'], 'preACZ_T1fromM0.nii.gz')
    subject['postACZ_T1fromM0_path'] = os.path.join(subject['ASLdir'], 'postACZ_T1fromM0.nii.gz')
    subject['postACZ_T1fromM0_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_T1fromM0_2preACZ.nii.gz')

    subject['preACZ_mask_path'] = os.path.join(subject['ASLdir'], 'preACZ_M0_brain_mask.nii.gz')
    subject['postACZ_mask_path'] = os.path.join(subject['ASLdir'], 'postACZ_M0_brain_mask.nii.gz')
    subject['postACZ_mask_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_M0_brain_mask_2preACZ.nii.gz')

    subject['preACZ_CBF_path'] = os.path.join(subject['ASLdir'], 'preACZ_QASL_2tolastPLD_forCBF/native_space/perfusion_calib.nii.gz')
    subject['postACZ_CBF_path'] = os.path.join(subject['ASLdir'], 'postACZ_QASL_2tolastPLD_forCBF/native_space/perfusion_calib.nii.gz')
    subject['postACZ_CBF_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_CBF_2preACZ.nii.gz')

    subject['preACZ_AAT_path'] = os.path.join(subject['ASLdir'], 'preACZ_QASL_allPLD_forAAT/native_space/arrival.nii.gz')
    subject['postACZ_AAT_path'] = os.path.join(subject['ASLdir'], 'postACZ_QASL_allPLD_forAAT/native_space/arrival.nii.gz')
    subject['postACZ_AAT_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_AAT_2preACZ.nii.gz')

    subject['preACZ_ATA_path'] = os.path.join(subject['ASLdir'], 'preACZ_QASL_1to2PLD_forATA/native_space/perfusion_calib.nii.gz')
    subject['postACZ_ATA_path'] = os.path.join(subject['ASLdir'], 'postACZ_QASL_1to2PLD_forATA/native_space/perfusion_calib.nii.gz')
    subject['postACZ_ATA_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_ATA_2preACZ.nii.gz')

    # Step 1: Convert DICOM to NIFTI, move input PACS DICOMS inputdir to SUBJECTdir/DICOMORIG/, and make this the new DICOMdir for further processing
    subject['DICOMdir'] = asl_convert_dicom_to_nifti(subject['DICOMdir'], subject['NIFTIdir'], imager='IMAGER', subject_dir=subject['SUBJECTdir'])

    # Step 2: Get SOURCE NIFTI files using string tags 
    filepre = sorted([f for f in os.listdir(subject['NIFTIdir']) if 'SOURCE_ASL' in f and 'preACZ' in f and f.endswith('2.nii.gz')])
    filepost = sorted([f for f in os.listdir(subject['NIFTIdir']) if 'SOURCE_ASL' in f and 'postACZ' in f and f.endswith('2.nii.gz')])
    if len(filepre) > 1 or len(filepost) > 1:
        print('WARNING: Multiple pre/post ASL datasets found. Using latest.')

    subject['preACZ_sourceNIFTI_path']  = filepre[-1]
    subject['postACZ_sourceNIFTI_path'] = filepost[-1]
    subject['preACZ_sourceDCM_path']    = subject['preACZ_sourceNIFTI_path'][:-7]
    subject['postACZ_sourceDCM_path']   = subject['postACZ_sourceNIFTI_path'][:-7]

    # Step 3: Locate template DICOMs
    templatedicom_files = os.listdir(subject['DICOMdir'])
    subject['preACZ_templateDCM_CBF_path'] = next(f for f in templatedicom_files if 'sWIP' in f and 'CBF' in f and 'preACZ' in f)
    subject['preACZ_templateDCM_AAT_path'] = next(f for f in templatedicom_files if 'sWIP' in f and 'AAT' in f and 'preACZ' in f)
    subject['preACZ_templateDCM_CVR_path'] = next(f for f in templatedicom_files if 'sWIP' in f and 'CVR' in f and 'preACZ' in f)
    subject['preACZ_templateDCM_ATA_path'] = next(f for f in templatedicom_files if 'sWIP' in f and 'ATA' in f and 'preACZ' in f)
   
    subject['postACZ_templateDCM_CBF_path'] = next(f for f in templatedicom_files if 'sWIP' in f and 'CBF' in f and 'postACZ' in f)
    subject['postACZ_templateDCM_AAT_path'] = next(f for f in templatedicom_files if 'sWIP' in f and 'AAT' in f and 'postACZ' in f)
    subject['postACZ_templateDCM_ATA_path'] = next(f for f in templatedicom_files if 'sWIP' in f and 'ATA' in f and 'postACZ' in f)

    # Step 4: DICOM info extraction and LL correction
    subject = asl_extract_params_dicom(subject, subject['preACZ_sourceDCM_path'])
    subject['preACZ_templateNII_path'] = os.path.join(subject['NIFTIdir'], subject['preACZ_sourceNIFTI_path']) 
    subject['postACZ_templateNII_path'] = os.path.join(subject['NIFTIdir'], subject['postACZ_sourceNIFTI_path']) 

    # Step 5: Compute Look Locker correction on Mxy per PLD   
    subject['LookLocker_correction_factor_perPLD'] = asl_look_locker_correction(subject)

    # Step 6: Split/control-label, save NIFTI
    subject = asl_prepare_asl_data(subject, subject['preACZ_sourceNIFTI_path'], 'preACZ')
    subject = asl_prepare_asl_data(subject, subject['postACZ_sourceNIFTI_path'], 'postACZ')

    # Step 7: Brain extraction for mask and T1 from M0
    subject = asl_bet_t1_from_m0(subject,'preACZ', 'fast') 
    subject = asl_bet_t1_from_m0(subject,'postACZ', 'fast')

    for prefix in ['preACZ', 'postACZ']:
    # Step 8: Motion Correction
        #asl_motioncorrection_ants(subject[f'{prefix}_PLDall_labelcontrol_path'], subject[f'{prefix}_M0_path'], subject[f'{prefix}_PLDall_labelcontrol_path'])
        #asl_motioncorrection_ants(subject[f'{prefix}_PLD2tolast_labelcontrol_path'], subject[f'{prefix}_M0_path'], subject[f'{prefix}_PLD2tolast_labelcontrol_path'])
        #asl_motioncorrection_ants(subject[f'{prefix}_PLD1to2_labelcontrol_path'], subject[f'{prefix}_M0_path'], subject[f'{prefix}_PLD1to2_labelcontrol_path'])

    # Step 9: ASL Quantification analysis
        # all PLD for AAT (arterial arrival time map)
        asl_qasl_analysis(subject, subject[f'{prefix}_PLDall_labelcontrol_path'], subject[f'{prefix}_M0_path'], subject[f'{prefix}_mask_path'] , os.path.join(subject['ASLdir'], f'{prefix}_QASL_allPLD_forAAT'), subject['PLDS'][0:], subject['inference_method'])
       
        # 2-to-last PLD for CBF map
        asl_qasl_analysis(subject, subject[f'{prefix}_PLD2tolast_labelcontrol_path'], subject[f'{prefix}_M0_path'], subject[f'{prefix}_mask_path'] , os.path.join(subject['ASLdir'], f'{prefix}_QASL_2tolastPLD_forCBF'), subject['PLDS'][1:], subject['inference_method'])
       
        # 1to2 PLDs for ATA map ->  then do no fit for the arterial component 'artoff'
        asl_qasl_analysis(subject, subject[f'{prefix}Z_PLD1to2_labelcontrol_path'], subject[f'{prefix}_M0_path'], subject[f'{prefix}_mask_path'] , os.path.join(subject['ASLdir'], f'{prefix}_QASL_1to2PLD_forATA'), subject['PLDS'][0:2], subject['inference_method'], artoff='artoff')
    
    # Step 10: register post-ACZ ASL data to pre-ACZ ASL data using Elastix 
    asl_registration_prepostACZ(subject)

    # Step 11: Generate CBF/AAT/ATA/CVR results (nifti, dicom PACS, .pngs) for pre- and postACZ, including registration of postACZ to preACZ as reference data, and target for computed CVR map
    subject = asl_save_results_cbfaatcvr(subject)

    print('-- Finished --')
    return subject
