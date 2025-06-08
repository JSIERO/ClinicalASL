
% ClinicalASL toolbox 2023, JCWSiero
%%%%%%%%%%%%%%%%%%%%% ASL Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by Jeroen Siero  25-05-2023 for the MOYAMOYA study
% includes automatic DICOM file loading, anatomy segmentation and  registration, outlier removal, data construction, BASIL analysis, CBF, map smoothing, CVR registration, calculation and saving
clear all
close all
clc

SUBJECT.GITHUB_ClinicalASLDIR = '/home/jeroen/GITHUB/ClinicalASL/';
SUBJECT.masterdir = '/Fridge/users/jeroen/MOYAMOYA/';
SUBJECT.RegistrationMethod = 'elastix'; %choose: 'matlab_imreg', or 'elastix'
addpath(SUBJECT.GITHUB_ClinicalASLDIR)
addpath([SUBJECT.GITHUB_ClinicalASLDIR 'MNI'])
addpath([SUBJECT.GITHUB_ClinicalASLDIR 'generalFunctions'])
addpath([SUBJECT.GITHUB_ClinicalASLDIR 'generalFunctions/export_fig'])
setenv('PATH', [getenv('PATH') ';' SUBJECT.GITHUB_ClinicalASLDIR 'generalFunctions'])

% global parameters
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
SUBJECT.range_AATdelta = [-0.7 0.7];% delta time (s), delta arterial arrival time, postACZ - vTR
SUBJECT.range_aCBV = [0 2]; % arterial blodo volume estimate in volume fraction voxel (%)

%% %%%%%%%%%%%%%%%%%%%%%%% 1. Subject information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get subject folder name, select folder containing all patient data

SUBJECT.SUBJECTdir = uigetdir(SUBJECT.masterdir,'Select subject folder');

% create folder paths
SUBJECT.DICOMdir = fullfile(SUBJECT.SUBJECTdir,'DICOM/'); % DICOM  path
SUBJECT.NIFTIdir = fullfile(SUBJECT.SUBJECTdir,'NIFTI/'); % NIFTI  path
SUBJECT.ASLdir = fullfile(SUBJECT.SUBJECTdir,'ASL_vTR/'); % ASL path
SUBJECT.RESULTSdir = fullfile(SUBJECT.SUBJECTdir,'ASL_vTR','FIGURE_RESULTS'); % RESULTS path

% create folders
if logical(max(~isfolder({SUBJECT.NIFTIdir; SUBJECT.ASLdir; SUBJECT.RESULTSdir})))
    mkdir(SUBJECT.NIFTIdir); % create NIFTI folder
    mkdir(SUBJECT.ASLdir); % create ASL folder
    mkdir(SUBJECT.RESULTSdir); % create RESULTS folder
end

% convert and rename DICOM files in DICOM folder to NIFTI folder
SUBJECT = ASLConvertDICOMtoNIFTI(SUBJECT,[], 'vTR_only')

% Get vTR-ASL nifti filenames
filenames_vTR = dir(fullfile(SUBJECT.NIFTIdir, '*PRIDE*SOURCE*vTR*.nii.gz'));% find SOURCE data ASL
filenamesize=zeros(1,length(filenames_vTR));
for i=1:length(filenames_vTR)
    filenamesize(i)=sum(filenames_vTR(i,1).name); % sum the string values to find fielname with smallest scannumber (=vTR)
end
[filenamesize_sort, Isort]=sort(filenamesize);
filevTR = filenames_vTR(Isort(1),1).name;
filepostACZ = filenames_vTR(Isort(2),1).name;
SUBJECT.vTRfilenameNIFTI = filevTR;
SUBJECT.vTRfilenameDCM = filevTR(1:end-7);
SUBJECT.dummyfilenameSaveNII = fullfile(SUBJECT.NIFTIdir, SUBJECT.vTRfilenameNIFTI); % location .nii.gz NIFTI to be used as dummy template for saving NII's in the tool

% extract scan parameters from DICOM
SUBJECT = ASLExtractParamsDICOM_vTR(SUBJECT, SUBJECT.vTRfilenameDCM);

% Get vTR-ASL M0 nifti filenames
filenames_vTRM0 = dir(fullfile(SUBJECT.NIFTIdir, '*SOURCE*M0*.nii.gz'));% find SOURCE data ASL
filenamesizeM0=zeros(1,length(filenames_vTRM0));
for i=1:length(filenames_vTRM0)
    filenamesizeM0(i)=sum(filenames_vTRM0(i,1).name); % sum the string values to find fielname with smallest scannumber (=vTR)
end
[filenamesizeM0_sort, Isort_M0]=sort(filenamesizeM0);
SUBJECT.vTRfilenameDCM_M0 = fullfile(SUBJECT.DICOMdir, filenames_vTRM0(Isort_M0(1),1).name(1:end-7)); % vTR M0 DICOM used to generate DICOMs of the analysis results in ASLSaveResultsCBFAATCVR_vTRASL_Windows()

SUBJECT.vTRM0filenameNIFTI = fullfile(SUBJECT.NIFTIdir, filenames_vTRM0(Isort_M0(1),1).name);

% CBF
SUBJECT.vTRCBFfilenameNIFTI = fullfile(SUBJECT.NIFTIdir, [SUBJECT.vTRfilenameNIFTI(1:end-7), '.nii.gz']);

% AAT
SUBJECT.vTRAATfilenameNIFTI = fullfile(SUBJECT.NIFTIdir, [SUBJECT.vTRfilenameNIFTI(1:end-7), 'a.nii.gz']);

SUBJECT.dummyfilenameSaveNII_M0 = SUBJECT.vTRM0filenameNIFTI; % location .nii.gz NIFTI to be used as dummy template for saving NII's in the tool
SUBJECT.dummyfilenameSaveNII_CBF = SUBJECT.vTRCBFfilenameNIFTI; % location .nii.gz NIFTI to be used as dummy template for saving NII's in the tool
SUBJECT.dummyfilenameSaveNII_AAT = SUBJECT.vTRAATfilenameNIFTI; % location .nii.gz NIFTI to be used as dummy template for saving NII's in the tool

disp('DICOMs converted to NIFTI');
%% %%%%%%%%%%%%%%%%%%%%%%%%  7. Generate resulting CBF/CVR/AAT/aCBV .png images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJECT = ASLSaveResultsCBFAAT_vTRASL(SUBJECT);

%% %%%%%%%%%%%%%%%%%%%%%%%% 8. Save workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Workspace_ClinicalASL_vTR = fullfile(SUBJECT.SUBJECTdir, 'Workspace_ClinicalASL_vTR.mat');
save(Workspace_ClinicalASL_vTR,'-v7.3');
disp('-- Finished --');

