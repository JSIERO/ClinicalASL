
% ClinicalASL toolbox 2023, JCWSiero
%%%%%%%%%%%%%%%%%%%%% ASL Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by Jeroen Siero  25-05-2023 for the MOYAMOYA study
% includes automatic DICOM file loading, anatomy segmentation and  registration, outlier removal, data construction, BASIL analysis, CBF, map smoothing, CVR registration, calculation and saving
clear all
close all
clc

SUBJECT.GITHUB_ClinicalASLDIR = 'J:\OneDrive\Documents\GitHub\ClinicalASL';
SUBJECT.masterdir = 'G:\DATA\MOYAMOYA\';
% addpath(SUBJECT.GITHUB_ClinicalASLDIR)
% addpath([SUBJECT.GITHUB_ClinicalASLDIR '\MNI'])
% addpath([SUBJECT.GITHUB_ClinicalASLDIR '\generalFunctions'])
% addpath([SUBJECT.GITHUB_ClinicalASLDIR '\generalFunctions\export_fig'])
setenv('PATH', [getenv('PATH') ';' SUBJECT.GITHUB_ClinicalASLDIR '\generalFunctions'])

% global parameters
SUBJECT.N_BS = 4; % Number of background suppression pulses
SUBJECT.labeleff = 0.85; %PCASL label efficiency
SUBJECT.lambda = 0.9; %water partition fraction, Alsop MRM 2014
SUBJECT.T1t = 1.3;    %s T1 of tissue DEFAULT is 1.3 at 3T, Alsop MRM 2014
SUBJECT.T1b = 1.65;   %s T1 of arterial blood DEFAULT is 1.65 at 3T, Alsop MRM 2014
SUBJECT.FWHM = 6; % smoothing kernel size 6mm FWHM, for CBF, AAT
SUBJECT.FWHM_M0 = 5; % smoothing kernel size  5mm FWHM, for M0_forQCBF for manual quantification
SUBJECT.range_adult_cbf = [0 75]; % colourbar range for adult CBF values
SUBJECT.range_child_cbf = [0 125]; % colourbar range for child CBF values
SUBJECT.range_cvr = [-50 50]; % colourbar range for CVR values
SUBJECT.range_AAT = [0.5 2.5]; % time (s), arterial arrival time
SUBJECT.range_AATdelta = [-1.2 1.2];% delta time (s), delta arterial arrival time, postACZ - preACZ
SUBJECT.range_aCBV = [0 2]; % arterial blodo volume estimate in volume fraction voxel (%)

%% %%%%%%%%%%%%%%%%%%%%%%% 1. Subject information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get subject folder name, select folder containing all patient data

SUBJECT.SUBJECTdir = uigetdir(SUBJECT.masterdir,'Select subject folder');

% create folder paths
SUBJECT.ANATOMYdir = [SUBJECT.SUBJECTdir,'\ANATOMY\']; % T1 anatomy path
SUBJECT.MNIdir = [SUBJECT.masterdir 'MNI\']; % MNI path, needs MNI_T1_2mm_brain MNI_BRAINMASK_2mm, and seg_0, seg_1, seg_2 (CSF, GM and WM) tissue segmentations: obtain from GITHUB\ClinicalASL
SUBJECT.SUBJECTMNIdir = [SUBJECT.SUBJECTdir '\MNI\']; % MNI path
SUBJECT.DICOMdir = [SUBJECT.SUBJECTdir,'\DICOM\']; % DICOM  path
SUBJECT.NIFTIdir = [SUBJECT.SUBJECTdir,'\NIFTI\']; % NIFTI  path
SUBJECT.ASLdir = [SUBJECT.SUBJECTdir,'\ASL_vTR\']; % ASL path
SUBJECT.RESULTSdir = [SUBJECT.SUBJECTdir,'\ASL_vTR\FIGURE_RESULTS\']; % RESULTS path

% create folders
if logical(max(~isfolder({SUBJECT.ANATOMYdir; SUBJECT.NIFTIdir; SUBJECT.ASLdir; SUBJECT.RESULTSdir; SUBJECT.SUBJECTMNIdir})))
    mkdir(SUBJECT.ANATOMYdir); % create Anatomy folder
    mkdir(SUBJECT.SUBJECTMNIdir); % create subject MNI folder
    mkdir(SUBJECT.NIFTIdir); % create NIFTI folder
    mkdir(SUBJECT.ASLdir); % create ASL folder
    mkdir(SUBJECT.RESULTSdir); % create RESULTS folder
end

% convert and rename DICOM files in DICOM folder to NIFTI folder
ASLConvertDICOMtoNIFTI(SUBJECT.DICOMdir, SUBJECT.NIFTIdir, 'vTR_only');

% Get vTR-ASL nifti filenames
filenames_vTR = dir([SUBJECT.NIFTIdir, '*PRIDE*SOURCE*vTR*.nii.gz']);% find SOURCE data ASL
filenamesize=zeros(1,length(filenames_vTR));
for i=1:length(filenames_vTR)
    filenamesize(i)=sum(filenames_vTR(i,1).name); % sum the string values to find fielname with smallest scannumber (=preACZ)
end
[filenamesize_sort, Isort]=sort(filenamesize);
filepreACZ = filenames_vTR(Isort(1),1).name;
filepostACZ = filenames_vTR(Isort(2),1).name;
SUBJECT.preACZfilenameNIFTI = filepreACZ;
SUBJECT.postACZfilenameNIFTI = filepostACZ;
SUBJECT.preACZfilenameDCM = filepreACZ(1:end-7);

% extract scan parameters from DICOM
SUBJECT = ASLExtractParamsDICOM_vTR(SUBJECT, SUBJECT.preACZfilenameDCM);

% Get vTR-ASL M0, CBF and AAT nifti filenames
filenames_vTRM0 = dir([SUBJECT.NIFTIdir, '*SOURCE*M0*.nii.gz']);% find SOURCE data ASL
filenamesizeM0=zeros(1,length(filenames_vTRM0));
for i=1:length(filenames_vTRM0)
    filenamesizeM0(i)=sum(filenames_vTRM0(i,1).name); % sum the string values to find fielname with smallest scannumber (=preACZ)
end
[filenamesizeM0_sort, Isort_M0]=sort(filenamesizeM0);
SUBJECT.preACZM0filenameNIFTI = [SUBJECT.NIFTIdir filenames_vTRM0(Isort_M0(1),1).name];
SUBJECT.postACZM0filenameNIFTI = [SUBJECT.NIFTIdir filenames_vTRM0(Isort_M0(2),1).name];

% CBF
SUBJECT.preACZCBFfilenameNIFTI = [SUBJECT.NIFTIdir SUBJECT.preACZfilenameNIFTI(1:end-7) '.nii.gz'];
SUBJECT.postACZCBFfilenameNIFTI = [SUBJECT.NIFTIdir SUBJECT.postACZfilenameNIFTI(1:end-7) '.nii.gz'];

% AAT
SUBJECT.preACZAATfilenameNIFTI = [SUBJECT.NIFTIdir SUBJECT.preACZfilenameNIFTI(1:end-7) 'a.nii.gz'];
SUBJECT.postACZAATfilenameNIFTI = [SUBJECT.NIFTIdir SUBJECT.postACZfilenameNIFTI(1:end-7) 'a.nii.gz'];

SUBJECT.dummyfilenameSaveNII_M0 = SUBJECT.preACZM0filenameNIFTI; % location .nii.gz NIFTI to be used as dummy template for saving NII's in the tool
SUBJECT.dummyfilenameSaveNII_CBF = SUBJECT.preACZCBFfilenameNIFTI; % location .nii.gz NIFTI to be used as dummy template for saving NII's in the tool
SUBJECT.dummyfilenameSaveNII_AAT = SUBJECT.preACZAATfilenameNIFTI; % location .nii.gz NIFTI to be used as dummy template for saving NII's in the tool

disp('DICOMs converted to NIFTI');

%% %%%%%%%%%%%%%%%%%%%%%%%%  7. Generate resulting CBF\CVR\AAT\aCBV .png images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJECT = ASLSaveResultsCBFAATCVR_vTRASL_Windows(SUBJECT);

%% %%%%%%%%%%%%%%%%%%%%%%%% 8. Save workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Workspace_ClinicalASL_vTR = [SUBJECT.SUBJECTdir '\Workspace_ClinicalASL_vTR.mat'];
save(Workspace_ClinicalASL_vTR,'-v7.3');
disp('-- Finished --');

