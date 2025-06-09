import os
from clinical_asl_pipeline.asl_convert_dicom_to_nifti import asl_convert_dicom_to_nifti
from clinical_asl_pipeline.asl_extract_params_dicom import  asl_extract_params_dicom
from clinical_asl_pipeline.asl_look_locker_correction import  asl_look_locker_correction
from clinical_asl_pipeline.asl_prepare_asl_data import  asl_prepare_asl_data
from clinical_asl_pipeline.asl_t1_from_m0_processing import  asl_t1_from_m0_processing
from clinical_asl_pipeline.asl_basil_analysis import  asl_basil_analysis
from clinical_asl_pipeline.asl_qasl_analysis import  asl_qasl_analysis
from clinical_asl_pipeline.asl_motioncorrection_ants import asl_motioncorrection_ants
from clinical_asl_pipeline.asl_save_results_cbfaatcvr import asl_save_results_cbfaatcvr

def mri_diamox_umcu_clinicalasl_cvr_imager(inputdir, outputdir):
    subject = {}
    subject['DICOMdir'] = inputdir # input folder for extracted PACS DICOM
    subject['SUBJECTdir'] = outputdir # inpput folder for results and generated DICOMS of ALS derived imaged for PACS
    docker_compiled_app_location = '/home/jeroen/GITHUB/ClinicalASL_Python/clinical_asl_pipeline/'
    
    # Load registration and BASIL settings
    subject['ElastixParameterFile'] = os.path.join(docker_compiled_app_location, 'Par0001rigid_6DOF_MI_NIFTIGZ.txt')
    if not os.path.isfile(subject['ElastixParameterFile']):
        print(f"Missing Elastix file: {subject['ElastixParameterFile']}")
        return

    subject['locationBASILinfo'] = os.path.join(docker_compiled_app_location, 'BASIL_OPTIONS.txt')
    if not os.path.isfile(subject['locationBASILinfo']):
        print(f"Missing BASIL_OPTIONS file: {subject['locationBASILinfo']}")
        return

    # Default parameters
    subject.update({
        'tau': 2,
        'N_BS': 4,
        'labeleff': 0.85,
        'lambda': 0.9,
        'T1t': 1.3,
        'T1b': 1.65,
        'FWHM': 6,
        'FWHM_M0': 5,
        'range_adult_cbf': [0, 75],
        'range_child_cbf': [0, 125],
        'range_cvr': [-50, 50],
        'range_AAT': [0.5, 2.5],
        'range_ATA': [0, 125]
    })

    # Paths
    subject['NIFTIdir'] = os.path.join(subject['SUBJECTdir'], 'NIFTI')
    subject['ASLdir'] = os.path.join(subject['SUBJECTdir'], 'ASL')
    subject['RESULTSdir'] = os.path.join(subject['ASLdir'], 'FIGURE_RESULTS')
    subject['DICOMRESULTSdir'] = subject['SUBJECTdir'] # here the output DICOMS with ALS derived images for PACS will be saved

    for path in [subject['NIFTIdir'], subject['ASLdir'], subject['RESULTSdir']]:
        os.makedirs(path, exist_ok=True)
        
    # Define all necessary file paths
    subject['preACZ_T1fromM0_path'] = os.path.join(subject['ASLdir'], 'preACZ_T1fromM0.nii.gz')
    subject['postACZ_T1fromM0_path'] = os.path.join(subject['ASLdir'], 'postACZ_T1fromM0.nii.gz')
    subject['postACZ_T1fromM0_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_T1fromM0_2preACZ.nii.gz')

    subject['preACZ_mask_path'] = os.path.join(subject['ASLdir'], 'preACZ_M0_brain_mask.nii.gz')
    subject['postACZ_mask_path'] = os.path.join(subject['ASLdir'], 'postACZ_M0_brain_mask.nii.gz')
    subject['postACZ_mask_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_M0_brain_mask_2preACZ.nii.gz')

    subject['preACZ_CBF_path'] = os.path.join(subject['ASLdir'], 'preACZ_BASIL_2tolastPLD_forCBF/native_space/perfusion_calib.nii.gz')
    subject['postACZ_CBF_path'] = os.path.join(subject['ASLdir'], 'postACZ_BASIL_2tolastPLD_forCBF/native_space/perfusion_calib.nii.gz')
    subject['postACZ_CBF_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_CBF_2preACZ.nii.gz')

    subject['preACZ_AAT_path'] = os.path.join(subject['ASLdir'], 'preACZ_BASIL_allPLD_forAAT/native_space/arrival.nii.gz')
    subject['postACZ_AAT_path'] = os.path.join(subject['ASLdir'], 'postACZ_BASIL_allPLD_forAAT/native_space/arrival.nii.gz')
    subject['postACZ_AAT_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_AAT_2preACZ.nii.gz')

    subject['preACZ_ATA_path'] = os.path.join(subject['ASLdir'], 'preACZ_BASIL_1to2PLD_forATA/native_space/perfusion_calib.nii.gz')
    subject['postACZ_ATA_path'] = os.path.join(subject['ASLdir'], 'postACZ_BASIL_1to2PLD_forATA/native_space/perfusion_calib.nii.gz')
    subject['postACZ_ATA_2preACZ_path'] = os.path.join(subject['ASLdir'], 'postACZ_ATA_2preACZ.nii.gz')

    # Step 1: Convert DICOM to NIFTI, move input PACS DICOMS in input dir in SUBJECTdir/DICOMORIG/, and make this the new DICOMdir for further processing
    subject['DICOMdir'] = asl_convert_dicom_to_nifti(subject['DICOMdir'], subject['NIFTIdir'], imager='IMAGER', subject_dir=subject['SUBJECTdir'])

    # Step 2: Get NIFTI files using string tags from the 
    filepre = sorted([f for f in os.listdir(subject['NIFTIdir']) if 'SOURCE_ASL' in f and 'preACZ' in f and f.endswith('2.nii.gz')])
    filepost = sorted([f for f in os.listdir(subject['NIFTIdir']) if 'SOURCE_ASL' in f and 'postACZ' in f and f.endswith('2.nii.gz')])
    if len(filepre) > 1 or len(filepost) > 1:
        print("WARNING: Multiple pre/post ASL datasets found. Using latest.")

    subject['preACZfilenameNIFTI'] = filepre[-1]
    subject['postACZfilenameNIFTI'] = filepost[-1]
    subject['preACZfilenameDCM'] = subject['preACZfilenameNIFTI'][:-7]
    subject['postACZfilenameDCM'] = subject['postACZfilenameNIFTI'][:-7]

    # Step 3: Locate dummy DICOMs
    dummy_files = os.listdir(subject['DICOMdir'])
    subject['preACZfilenameDCM_CBF'] = next(f for f in dummy_files if 'sWIP' in f and 'CBF' in f and 'preACZ' in f)
    subject['preACZfilenameDCM_AAT'] = next(f for f in dummy_files if 'sWIP' in f and 'AAT' in f and 'preACZ' in f)
    subject['preACZfilenameDCM_CVR'] = next(f for f in dummy_files if 'sWIP' in f and 'CVR' in f and 'preACZ' in f)
    subject['preACZfilenameDCM_ATA'] = next(f for f in dummy_files if 'sWIP' in f and 'ATA' in f and 'preACZ' in f)
   
    subject['postACZfilenameDCM_CBF'] = next(f for f in dummy_files if 'sWIP' in f and 'CBF' in f and 'postACZ' in f)
    subject['postACZfilenameDCM_AAT'] = next(f for f in dummy_files if 'sWIP' in f and 'AAT' in f and 'postACZ' in f)
    subject['postACZfilenameDCM_ATA'] = next(f for f in dummy_files if 'sWIP' in f and 'ATA' in f and 'postACZ' in f)

    # Step 4: DICOM info extraction and LL correction
    subject = asl_extract_params_dicom(subject, subject['preACZfilenameDCM'])
    subject['dummyfilenameSaveNII'] = os.path.join(subject['NIFTIdir'], subject['preACZfilenameNIFTI']) 

    # Step 5: Compute Look Locker correction on Mxy per PLD   
    subject['LookLocker_correction_factor_perPLD'] = asl_look_locker_correction(subject)

    # Step 6: Split/control-label, save NIFTI
    subject = asl_prepare_asl_data(subject, subject['preACZfilenameNIFTI'], 'preACZ', 'fast')
    subject = asl_prepare_asl_data(subject, subject['postACZfilenameNIFTI'], 'postACZ', 'fast')
    print("DICOMs converted to NIFTI")

    # Step 7: Create T1 from M0
    print("Create T1 from M0")    
    subject = asl_t1_from_m0_processing(subject,'preACZ', 'fast') 
    subject = asl_t1_from_m0_processing(subject,'postACZ', 'fast')

    print("Perform ASL Quantification analysis")
    for prefix in ['preACZ', 'postACZ']:
        m0 = os.path.join(subject['ASLdir'], f'{prefix}_M0.nii.gz')
        mask = os.path.join(subject['ASLdir'], f'{prefix}_M0_brain_mask.nii.gz')    
    # Step 8: Motion Correction
        print("Perform motion correction")
        PLDall = os.path.join(subject['ASLdir'], f'{prefix}_allPLD_label1label2.nii.gz') 
        PLD2tolast = os.path.join(subject['ASLdir'], f'{prefix}_2tolastPLD_label1label2.nii.gz')
        PLD1to2 = os.path.join(subject['ASLdir'], f'{prefix}_1to2PLD_label1label2.nii.gz')
        asl_motioncorrection_ants(PLDall, m0, PLDall)
        asl_motioncorrection_ants(PLD2tolast, m0, PLD2tolast)
        asl_motioncorrection_ants(PLD1to2, m0, PLD1to2)

    # Step 9: ASL Quantification analysis
    # all PLD for AAT (arterial arrival time map)
        #asl_basil_analysis(subject, PLDall, m0, mask, os.path.join(subject['ASLdir'], f'{prefix}_BASIL_allPLD_forAAT'), subject['PLDS'][0:], subject['locationBASILinfo'])
        asl_qasl_analysis(subject, PLDall, m0, mask, os.path.join(subject['ASLdir'], f'{prefix}_BASIL_allPLD_forAAT'), subject['PLDS'][0:], 'ssvb')
       
    # 2-to-last PLD for CBF map
        #asl_basil_analysis(subject, PLD2tolast, m0, mask, os.path.join(subject['ASLdir'], f'{prefix}_BASIL_2tolastPLD_forCBF'), subject['PLDS'][1:], subject['locationBASILinfo'])
        asl_qasl_analysis(subject, PLD2tolast, m0, mask, os.path.join(subject['ASLdir'], f'{prefix}_BASIL_2tolastPLD_forCBF'), subject['PLDS'][1:], 'ssvb')
       
    # 1to2 PLDs for ATA map ->  then do no fit for the arterial component 'artoff'
        #asl_basil_analysis(subject, PLD1to2, m0, mask, os.path.join(subject['ASLdir'], f'{prefix}_BASIL_1to2PLD_forATA'),subject['PLDS'][0:2],subject['locationBASILinfo'], 'artoff')
        asl_qasl_analysis(subject, PLD1to2, m0, mask, os.path.join(subject['ASLdir'], f'{prefix}_BASIL_1to2PLD_forATA'), subject['PLDS'][0:2], 'ssvb', artoff='artoff')
    
    # Step 9: Generate CBF/AAT/ATA/CVR results (nifti, dicom PACS, .pngs) for pre- and postACZ, including registration of postACZ to preACZ as reference data, and target for computed CVR map
    subject = asl_save_results_cbfaatcvr(subject)

    print("-- Finished --")
    return subject
