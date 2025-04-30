
% ClinicalASL toolbox 2023, JCWSiero
%%%%%%%%%%%%%%%%%%%%% ASL Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by Jeroen Siero  25-05-2023 for the MOYAMOYA study
% includes automatic DICOM file loading, anatomy segmentation and  registration, outlier removal, data construction, BASIL analysis, CBF, map smoothing, CVR registration, calculation and saving
clear all
close all
clc
SUBJECT.masterdir='/Fridge/users/jeroen/MOYAMOYA/';

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

%% %%%%%%%%%%%%%%%%%%%%%%% 1. Subject information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get subject folder name, select folder containing all patient data

SUBJECT.SUBJECTdir = uigetdir(SUBJECT.masterdir,'Select subject folder');

% create folder paths
SUBJECT.DICOMdir = [SUBJECT.SUBJECTdir,'/DICOM/']; % DICOM  path
SUBJECT.NIFTIdir = [SUBJECT.SUBJECTdir,'/NIFTI/']; % NIFTI  path
SUBJECT.ASLdir = [SUBJECT.SUBJECTdir,'/ASL/']; % ASL path
SUBJECT.RESULTSdir = [SUBJECT.SUBJECTdir,'/ASL/FIGURE_RESULTS/']; % RESULTS path
% extra FSL BASIL options .txt location
SUBJECT.locationBASILinfo=[SUBJECT.masterdir 'BASIL_OPTIONS.txt']; % location .txt file with addition model options for CBF quantification BASIL
if ~isfile(SUBJECT.locationBASILinfo)
    warning('no BASIL_OPTIONS.txt file found in study folder (masterdir), please copy from GITHUB/ClinicalASL/ repository')
    return
end
% create folders
if logical(max(~isfolder({SUBJECT.NIFTIdir; SUBJECT.ASLdir; SUBJECT.RESULTSdir})))
    mkdir(SUBJECT.NIFTIdir); % create NIFTI folder
    mkdir(SUBJECT.ASLdir); % create ASL folder
    mkdir(SUBJECT.RESULTSdir); % create RESULTS folder
end

% convert and rename DICOM files in DICOM folder to NIFTI folder
ASLConvertDICOMtoNIFTI(SUBJECT.DICOMdir, SUBJECT.NIFTIdir)

% Get ASL nifti filenames
% preACZ path
filepreACZ = dir([SUBJECT.NIFTIdir, '*SOURCE_ASL*preACZ*2.nii.gz']);% find SOURCE data ASL
% postACZ path
filepostACZ = dir([SUBJECT.NIFTIdir, '*SOURCE_ASL*postACZ*2.nii.gz']);% find SOURCE data ASL

if (size(filepreACZ,1) > 1 ) || (size(filepostACZ,1) > 1)
    warning(' More than 1 preAZC or postACZ ASL dataset found -  !! taking the last scanned dataset... !!')
end

SUBJECT.preACZfilenameNIFTI = filepreACZ(end,1).name;
SUBJECT.preACZfilenameDCM = filepreACZ(end,1).name(1:end-7);
SUBJECT.postACZfilenameNIFTI = filepostACZ(end,1).name;
SUBJECT.postACZfilenameDCM = filepostACZ(end,1).name(1:end-7);

%% %%%%%%%%%%%%%%%%%%%%%%% 2. Extract DICOM information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch scan parameters
SUBJECT = ASLExtractParamsDICOM(SUBJECT, SUBJECT.preACZfilenameDCM);
SUBJECT.dummyfilenameSaveNII = [SUBJECT.NIFTIdir SUBJECT.preACZfilenameNIFTI]; % lication .nii.gz NIFTI to be used as dummy template for saving NII's in the tool
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

%% %%%%%%%%%%%%%%%%%%%%%%%% 6. BASIL CBF Analysis for both Original and Outlier removed ASL data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Perform BASIL analysis for both original and outlier removed ASL data')

session = {'preACZ', 'postACZ'};

for i=1:length(session)
    prefix = char(session(i));
    % %%%% % all PLD for AAT (arterial arrival time map) % %%%%%
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_allPLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_allPLD_forAAT'], [1:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
    % %%%% % 2tolast PLD for CBF map % %%%% %
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_2tolastPLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_2tolastPLD_forCBF'], [2:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%  7. Generate resulting CBF/CVR/AAT/aCBV .png images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJECT = ASLSaveResultsCBFAATCVR_FAST(SUBJECT); %

%% %%%%%%%%%%%%%%%%%%%%%%%% 8. Save workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Workspace_ClinicalASL = [SUBJECT.SUBJECTdir '/Workspace_ClinicalASL.mat'];
save(Workspace_ClinicalASL,'-v7.3');
disp('-- Finished --');

