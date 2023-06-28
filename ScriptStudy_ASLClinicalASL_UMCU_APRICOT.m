%%%%%%%%%%%%%%%%%%%%% ASL Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% written by Jeroen Siero  25-05-2023 for the APRICOT study
%%% includes automatic DICOM file loading, anatomy segmentation and  registration, outlier removal, data construction, BASIL analysis, CBF, map smoothing,, calculation and saving
clear all
close all
clc
%% SOURCE DATA SUBJECTS
SUBJECT.masterdir='/Fridge/users/jeroen/APRICOT/';

SUBJECT.tau = 1.65; % Label duration
SUBJECT.N_BS = 2; % Number of background suppression pulses
SUBJECT.labeleff = 0.85; %PCASL label efficiency
SUBJECT.lambda = 0.9;%water partition fraction, Alsop MRM 2014
SUBJECT.T1t = 1.3;    %s T1 of tissue DEFAULT is 1.3 at 3T, Alsop MRM 2014
SUBJECT.T1b = 1.65;   %s T1 of arterial blood DEFAULT is 1.65 at 3T, Alsop MRM 2014
SUBJECT.FWHM = 4; % smoothing kernel size 6mm FWHM, for CBF, AAT
SUBJECT.FWHM_M0 = 5; % smoothing kernel size  5mm FWHM, for M0_forQCBF for manual quantification
SUBJECT.outlierFactor = 2.5; % outlierFactor x stds from mean CBF are considered outliers in step1 of Dolui et al SCORE method
SUBJECT.range_cbf = [0 100]; % colourbar range for CBF values
SUBJECT.range_cvr = [-50 50]; % colourbar range for CVR values
SUBJECT.range_AAT = [0.5 2.5]; % time (s), arterial arrival time
SUBJECT.range_aCBV = [0 2]; % arterial blodo volume estimate in volume fraction voxel (%)

% loop over subjects

subjnames_baseline=importdata([SUBJECT.masterdir 'subjectlist_baseline']);
subjnames_followup=importdata([SUBJECT.masterdir 'subjectlist_followup']);
subjnames=[subjnames_baseline; subjnames_followup];

subjnames={'APRICPOT_MR22_16chHEAD_17012023_ClinicalASLtest'}

%% %%%%%%%%%%%%%%%%%%%%%%% 1. Subject information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get subject folder name, select folder containing all patient data

SUBJECT.SUBJECTdir = [SUBJECT.masterdir char(subjnames)];

% create folder paths
SUBJECT.ANATOMYdir = [SUBJECT.SUBJECTdir,'/ANATOMY/']; % T1 anatomy path
SUBJECT.MNIdir = [SUBJECT.masterdir 'MNI/']; % MNI path, needs MNI_T1_2mm_brain MNI_BRAINMASK_2mm, and seg_0, seg_1, seg_2 (CSF, GM and WM) tissue segmentations
SUBJECT.SUBJECTMNIdir = [SUBJECT.SUBJECTdir '/MNI/']; % MNI path
SUBJECT.DICOMdir = [SUBJECT.SUBJECTdir,'/DICOM/']; % DICOM  path
SUBJECT.NIFTIdir = [SUBJECT.SUBJECTdir,'/NIFTI/']; % NIFTI  path
SUBJECT.ASLdir = [SUBJECT.SUBJECTdir,'/ASL/']; % ASL path
SUBJECT.RESULTSdir = [SUBJECT.SUBJECTdir,'/RESULTS/']; % RESULTS path
% extra FSL BASIL options .txt location
SUBJECT.locationBASILinfo=[SUBJECT.masterdir 'BASIL_OPTIONS.txt']; % location .txt file with addition model options for CBF quantification BASIL

% create folders
if logical(max(~isfolder({SUBJECT.ANATOMYdir; SUBJECT.DICOMdir; SUBJECT.NIFTIdir; SUBJECT.ASLdir; SUBJECT.RESULTSdir})))
    mkdir(SUBJECT.ANATOMYdir); % create Anatomy folder
    mkdir(SUBJECT.SUBJECTMNIdir); % create subject MNI folder
    mkdir(SUBJECT.DICOMdir); % create DICOM folder
    mkdir(SUBJECT.NIFTIdir); % create NIFTI folder
    mkdir(SUBJECT.ASLdir); % create ASL folder
    mkdir(SUBJECT.RESULTSdir); % create RESULTS folder
end

% convert and rename DICOM files in DICOM folder to NIFTI folder
ASLConvertDICOMtoNIFTI(SUBJECT.DICOMdir, SUBJECT.NIFTIdir)

% Get ASL nifti filenames
% SOURCE ASL NIFTI path
fileSOURCEASLNIFTI = dir([SUBJECT.NIFTIdir, '*SOURCE_ASL*2.nii.gz']); % find SOURCE data ASL

if size(fileSOURCEASLNIFTI,1) > 1
    warning(' More than 1 SOURCE ASL dataset found -  !! taking the last scanned dataset... !!')
end

SUBJECT.filenameNIFTI = fileSOURCEASLNIFTI(end,1).name;
SUBJECT.filenameDCM = fileSOURCEASLNIFTI(end,1).name(1:end-7);

%% %%%%%%%%%%%%%%%%%%%%%%%% 2. Extract DICOM information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch scan parameters
SUBJECT = ASLExtractParamsDICOM(SUBJECT, SUBJECT.filenameDCM);
SUBJECT.dummyfilenameSaveNII = [SUBJECT.NIFTIdir SUBJECT.filenameNIFTI]; % lication .nii.gz NIFTI to be used as dummy template for saving NII's in the tool
   
% Obtain Look-Locker correction factor
SUBJECT.LookLocker_correction_factor_perPLD = ASLLookLockerCorrectionFactor_mDelayPCASL(SUBJECT); % LookLocker correction factor, depending on the flipangle and PLDs

%% %%%%%%%%%%%%%%%%%%%%%%%% 3. Modify NIFTI to correct names, correct Mz loss (small fip angle) using Look Locker correction %%%%%%%%%%%%%%%%%%%%%
% save per PLD and control and label volumes (interleaved), and save all ASL and M0 in struct SUBJECT
SUBJECT = ASLPrepareASLData(SUBJECT, SUBJECT.filenameNIFTI, 'ASL');

%% %%%%%%%%%%%%%%%%%%%%%%%% 4. Generate T1 from M0 , T1 Tissue segmentation and registration to T1 anatomy and MNI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Locate T1 anatomy NIFTI
SUBJECT.T1ANATfilenameNIFTI = dir([SUBJECT.NIFTIdir, '*T1*3D*TFE*.nii*']); % find T1 anatomy NIFTI filename
SUBJECT.T1ANATfilenameNIFTI = SUBJECT.T1ANATfilenameNIFTI.name;

% T1w ANATOMY scan brain extraction and FSL FAST segmentation in GM, WM and CSF of  into ANATOMY dir
T1Processing(SUBJECT, SUBJECT.T1ANATfilenameNIFTI);

% create T1fromM0 and save M0 from ASL multiPLD data and tissue segmentation using FSL FAST
SUBJECT = ASLT1fromM0Processing(SUBJECT, 'ASL');

% Register T1 and tissue segmentations to ASL space
ASLT1Registration(SUBJECT,'ASL');

% Register MNI segmentations to ASL space using the T1fromM0
ASLMNIRegistration(SUBJECT,'ASL');
disp('T1 ASL MNI tissue segmentation and registration finished');

%% %%%%%%%%%%%%%%%%%%%%%%%% 5. Outlier identification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dolui et al. SCORE outlier method
tic
SUBJECT = ASLOutlierRemoval(SUBJECT, 'ASL');
toc
%% %%%%%%%%%%%%%%%%%%%%%%%% 6. BASIL CBF Analysis for both Original and Outlier removed ASL data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Perform BASIL analysis for both original and outlier removed ASL data')

prefix = 'ASL';
% %%%% % all PLD for AAT (arterial arrival time map) % %%%%%
ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_allPLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_allPLD_forAAT'], [1:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_allPLD_label1label2_OR.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_allPLD_forAAT_OR'], [1:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
% %%%% % 2tolast PLD for CBF map % %%%% %
ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_2tolastPLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_2tolastPLD_forCBF'], [2:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_2tolastPLD_label1label2_OR.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_2tolastPLD_forCBF_OR'], [2:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
% %%%% % 1to2 PLDs for ATA map -> compute COV --> then do no fit for the arterial component 'artoff' % %%%% %
ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_1to2PLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_1to2PLD_forATA'], [1:2], SUBJECT.locationBASILinfo,'artoff')
ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_1to2PLD_label1label2_OR.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_1to2PLD_forATA_OR'], [1:2], SUBJECT.locationBASILinfo,'artoff')
% %%%% % 1to2 PLDs for aCBV map  then do fit for the arterial component 'noartoff'%  %%%% %
ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_1to2PLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_1to2PLD_foraCBV'], [1:2], SUBJECT.locationBASILinfo, 'noartoff')
ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_1to2PLD_label1label2_OR.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_1to2PLD_foraCBV_OR'], [1:2], SUBJECT.locationBASILinfo, 'noartoff')

%% %%%%%%%%%%%%%%%%%%%%%%%%% 7. Generate resulting CBF/CVR/AAT/aCBV .png images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJECT = ASLSaveResultsCBFAAT(SUBJECT,'ASL', '');
SUBJECT = ASLSaveResultsCBFAAT(SUBJECT,'ASL', '_OR'); % save outlier removed data

%% %%%%%%%%%%%%%%%%%%%%%%%% 8. Save workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Workspace_ClinicalASL = [SUBJECT.SUBJECTdir '/Workspace_ClinicalASL.mat'];
save(Workspace_ClinicalASL,'-v7.3');
disp('-- Finished --');

