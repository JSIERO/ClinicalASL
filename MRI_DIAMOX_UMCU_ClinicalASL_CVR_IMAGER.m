function MRI_DIAMOX_UMCU_ClinicalASL_CVR_IMAGER(inputdir, outputdir)
% ClinicalASL toolbox 2025, JCWSiero
% for MRI DIAMOX scans, via IMAGER
% includes automatic DICOM file loading, anatomy segmentation and  registration, outlier removal, data construction, BASIL analysis, CBF, map smoothing, CVR registration, calculation and saving

%% %%%%%%%%%%%%%%%%%%%%%%% 1. Subject information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get subject folder name, select folder containing all patient data
SUBJECT.DICOMdir = inputdir; % input folder for extracted PACS DICOM
SUBJECT.SUBJECTdir = outputdir; 
docker_compiled_app_location = '/app/compiled_matlab_app/';

% location of registration Elastix File and FSL BASIL options
SUBJECT.ElastixParameterFile = fullfile(docker_compiled_app_location,'Par0001rigid_6DOF_MI_NIFTIGZ.txt'); % use 6DOF, rigidbody, Mutual information for registration
if ~isfile(SUBJECT.ElastixParameterFile)
    warning(['no Elastix registration TXT file found in ' SUBJECT.ElastixParameterFile ' , please copy from GITHUB/ClinicalASL/ repository'])
    return
end

SUBJECT.locationBASILinfo = fullfile(docker_compiled_app_location,'BASIL_OPTIONS.txt'); % location .txt file with addition model options for CBF quantification FSL BASIL
if ~isfile(SUBJECT.locationBASILinfo)
    warning(['no BASIL_OPTIONS.txt file found in ' SUBJECT.locationBASILinfo ', please copy from GITHUB/ClinicalASL/ repository'])
    return
end

% set parameters
SUBJECT.RegistrationMethod = 'elastix'; %choose: 'matlab_imreg', or 'elastix'
SUBJECT.tau = 2; % Label duration
SUBJECT.N_BS = 4; % Number of background suppression pulses
SUBJECT.labeleff = 0.85; %PCASL label efficiency
SUBJECT.lambda = 0.9;%water partition fraction, Alsop MRM 2014
SUBJECT.T1t = 1.3;    %s T1 of tissue DEFAULT is 1.3 at 3T, Alsop MRM 2014
SUBJECT.T1b = 1.65;   %s T1 of arterial blood DEFAULT is 1.65 at 3T, Alsop MRM 2014
SUBJECT.FWHM = 6; % smoothing kernel size 6mm FWHM, for CBF, AAT
SUBJECT.FWHM_M0 = 5; % smoothing kernel size  5mm FWHM, for M0_forQCBF for manual quantification
SUBJECT.range_adult_cbf = [0 75]; % colourbar range for adult CBF values
SUBJECT.range_child_cbf = [0 125]; % colourbar range for child CBF values
SUBJECT.range_cvr = [-50 50]; % colourbar range for CVR values
SUBJECT.range_AAT = [0.5 2.5]; % time (s), arterial arrival time

% create folder paths
SUBJECT.NIFTIdir = fullfile(SUBJECT.SUBJECTdir,'/NIFTI/'); % NIFTI  path
SUBJECT.ASLdir = fullfile(SUBJECT.SUBJECTdir,'/ASL/'); % ASL path
SUBJECT.RESULTSdir = fullfile(SUBJECT.SUBJECTdir,'/ASL/FIGURE_RESULTS/'); % RESULTS path
SUBJECT.DICOMRESULTSdir = fullfile(SUBJECT.SUBJECTdir); % DICOM RESULTS path

% set paths to files for registration (preACZ, post ACZ) and output
SUBJECT.preACZ_T1fromM0_path = fullfile(SUBJECT.ASLdir,'preACZ_T1fromM0.nii.gz');
SUBJECT.postACZ_T1fromM0_path = fullfile(SUBJECT.ASLdir, 'postACZ_T1fromM0.nii.gz');
SUBJECT.postACZ_T1fromM0_2preACZ_path = fullfile(SUBJECT.ASLdir, 'postACZ_T1fromM0_2preACZ.nii.gz');

SUBJECT.preACZ_mask_path = fullfile(SUBJECT.ASLdir, 'preACZ_M0_brain_mask.nii.gz');
SUBJECT.postACZ_mask_path = fullfile(SUBJECT.ASLdir, 'postACZ_M0_brain_mask.nii.gz');
SUBJECT.postACZ_mask_2preACZ_path = fullfile(SUBJECT.ASLdir, 'postACZ_M0_brain_mask_2preACZ.nii.gz');

SUBJECT.preACZ_CBF_path = fullfile(SUBJECT.ASLdir, 'preACZ_BASIL_2tolastPLD_forCBF', '/native_space/perfusion_calib.nii.gz');
SUBJECT.postACZ_CBF_path = fullfile(SUBJECT.ASLdir, 'postACZ_BASIL_2tolastPLD_forCBF', '/native_space/perfusion_calib.nii.gz');
SUBJECT.postACZ_CBF_2preACZ_path = fullfile(SUBJECT.ASLdir, 'postACZ_CBF_2preACZ.nii.gz');

SUBJECT.preACZ_AAT_path = fullfile(SUBJECT.ASLdir, 'preACZ_BASIL_allPLD_forAAT', '/native_space/arrival.nii.gz');
SUBJECT.postACZ_AAT_path = fullfile(SUBJECT.ASLdir, 'postACZ_BASIL_allPLD_forAAT', '/native_space/arrival.nii.gz');
SUBJECT.postACZ_AAT_2preACZ_path = fullfile(SUBJECT.ASLdir, 'postACZ_AAT_2preACZ.nii.gz');

SUBJECT.preACZ_ATA_path = fullfile(SUBJECT.ASLdir, 'preACZ_BASIL_1to2PLD_forATA', '/native_space/perfusion_calib.nii.gz');
SUBJECT.postACZ_ATA_path = fullfile(SUBJECT.ASLdir, 'postACZ_BASIL_1to2PLD_forATA', '/native_space/perfusion_calib.nii.gz');
SUBJECT.postACZ_ATA_2preACZ_path = fullfile(SUBJECT.ASLdir, 'postACZ_ATA_2preACZ.nii.gz');

% create folders
if logical(max(~isfolder({SUBJECT.NIFTIdir; SUBJECT.ASLdir; SUBJECT.RESULTSdir})))
    mkdir(SUBJECT.NIFTIdir); % create NIFTI folder
    mkdir(SUBJECT.ASLdir); % create ASL folder
    mkdir(SUBJECT.RESULTSdir); % create RESULTS folder
end

% convert and rename DICOM files in DICOM folder to NIFTI folder
SUBJECT = ASLConvertDICOMtoNIFTI(SUBJECT,'IMAGER')

% Get ASL nifti filenames

% preACZ path
filepreACZ = dir([SUBJECT.NIFTIdir, '*SOURCE_ASL*preACZ*2.nii.gz']);% find SOURCE data ASL
% postACZ path
filepostACZ = dir([SUBJECT.NIFTIdir, '*SOURCE_ASL*postACZ*2.nii.gz']);% find SOURCE data ASL

if (size(filepreACZ,1) > 1 ) || (size(filepostACZ,1) > 1)
    warning(' More than 1 preAZC or postACZ ASL dataset found -  !! taking the last scanned dataset... !!')
end

SUBJECT.preACZfilenameNIFTI = filepreACZ(end,1).name;
SUBJECT.postACZfilenameNIFTI = filepostACZ(end,1).name;
% Get ASL DICOM filenames
SUBJECT.preACZfilenameDCM = filepreACZ(end,1).name(1:end-7); % DICOM ASL source file
SUBJECT.postACZfilenameDCM = filepostACZ(end,1).name(1:end-7); % DICOM ASL source file

preACZfilenameDCM_CBF = dir(fullfile(SUBJECT.DICOMdir, 'sWIP*CBF*preACZ*'));% find dummy data CBF preACZ
preACZfilenameDCM_AAT = dir(fullfile(SUBJECT.DICOMdir, 'sWIP*AAT*preACZ*'));% find dummy data AAT preACZ
preACZfilenameDCM_CVR = dir(fullfile(SUBJECT.DICOMdir, 'sWIP*CVR*preACZ*'));% find dummy data CVR preACZ
preACZfilenameDCM_ATA = dir(fullfile(SUBJECT.DICOMdir, 'sWIP*ATA*preACZ*'));% find dummy data ATA preACZ

postACZfilenameDCM_CBF = dir(fullfile(SUBJECT.DICOMdir, 'sWIP*CBF*postACZ*'));% find dummy data CBF postACZ
postACZfilenameDCM_AAT = dir(fullfile(SUBJECT.DICOMdir, 'sWIP*AAT*postACZ*'));% find dummy data AAT postACZ
postACZfilenameDCM_ATA = dir(fullfile(SUBJECT.DICOMdir, 'sWIP*ATA*postACZ*'));% find dummy data ATA postACZ

SUBJECT.preACZfilenameDCM_CBF = preACZfilenameDCM_CBF.name; % DICOM CBF dummy file preACZ
SUBJECT.preACZfilenameDCM_AAT = preACZfilenameDCM_AAT.name; % DICOM AAT dummy file preACZ
SUBJECT.preACZfilenameDCM_CVR = preACZfilenameDCM_CVR.name; % DICOM CVR dummy file preACZ
SUBJECT.preACZfilenameDCM_ATA = preACZfilenameDCM_ATA.name; % DICOM ATA dummy file preACZ

SUBJECT.postACZfilenameDCM_CBF = postACZfilenameDCM_CBF.name; % DICOM CBF dummy file postACZ
SUBJECT.postACZfilenameDCM_AAT = postACZfilenameDCM_AAT.name; % DICOM AAT dummy file postACZ
SUBJECT.postACZfilenameDCM_ATA = postACZfilenameDCM_ATA.name; % DICOM CVR dummy file postACZ


%% %%%%%%%%%%%%%%%%%%%%%%% 2. Extract DICOM information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch scan parameters
SUBJECT = ASLExtractParamsDICOM(SUBJECT, SUBJECT.preACZfilenameDCM);
SUBJECT.dummyfilenameSaveNII = fullfile(SUBJECT.NIFTIdir, SUBJECT.preACZfilenameNIFTI); % location .nii.gz NIFTI to be used as dummy template for saving NII's in the tool
% Obtain Look-Locker correction factor
SUBJECT.LookLocker_correction_factor_perPLD = ASLLookLockerCorrectionFactor_mDelayPCASL(SUBJECT); % LookLocker correction factor, depending on the flipangle and PLDs

%% %%%%%%%%%%%%%%%%%%%%%%% 3. Modify NIFTI to correct names, correct Mz loss (small fip angle) using Look Locker correction %%%%%%%%%%%%%%%%%%%%%
% save per PLD and control and label volumes (interleaved), and save all ASL and M0 in struct SUBJECT
SUBJECT = ASLPrepareASLDataDICOM(SUBJECT, SUBJECT.preACZfilenameNIFTI, 'preACZ', 'fast'); % preACZ
SUBJECT = ASLPrepareASLDataDICOM(SUBJECT, SUBJECT.postACZfilenameNIFTI, 'postACZ', 'fast'); % postACZ

disp('DICOMs converted to NIFTI');

%% %%%%%%%%%%%%%%%%%%%%%%%% 4. Generate T1 from M0 , T1 Tissue segmentation and registration to T1 anatomy and MNI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create T1fromM0 from ASL multiPLD data
SUBJECT = ASLT1fromM0Processing(SUBJECT, 'preACZ', 'fast');
SUBJECT = ASLT1fromM0Processing(SUBJECT, 'postACZ','fast');

%% %%%%%%%%%%%%%%%%%%%%%%%% 5. BASIL CBF Analysis for both Original and Outlier removed ASL data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Perform BASIL analysis for both original and outlier removed ASL data')

session = {'preACZ', 'postACZ'};

for i=1:length(session)
    prefix = char(session(i));
    % %%%% % all PLD for AAT (arterial arrival time map) % %%%%%
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_allPLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_allPLD_forAAT'], [1:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
    % %%%% % 2tolast PLD for CBF map % %%%% %
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_2tolastPLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_2tolastPLD_forCBF'], [2:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
    % %%%% % 1to2 PLDs for ATA map ->  then do no fit for the arterial component 'artoff' % %%%% %
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_1to2PLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_1to2PLD_forATA'], [1:2], SUBJECT.locationBASILinfo,'artoff')
end

%% %%%%%%%%%%%%%%%%%%%%%%%% 6. Generate resulting CBF/CVR/AAT/aCBV .png images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJECT = ASLSaveResultsCBFAATCVR_FAST(SUBJECT); %

disp('-- Finished --');

