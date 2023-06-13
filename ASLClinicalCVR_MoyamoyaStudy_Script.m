%%%%%%%%%%%%%%%%%%%%% ASL Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% written by Jeroen Siero  25-05-2023
%%% modified by Jeroen Siero 26-11-2020, to yield correct small flip angle correction and M0 for the Look Lock multiPLD Moyamoya ASL scan
%%% modified by Jeroen Siero 09-12-2020, 2 compartement fit including macrovascular compartment
%%% modified by Jeroen Siero 04-01-2021, correct Mz due to loss caused by Look Locker readout, fixed tons of little things
%%% modified by Jeroen Siero 04-04-2023, complete rehaul of all scripts for efficiency, used now DICOM input 
%%% includes automatic DICOM file loading, anatomy segmentation and  registration, outlier removal, data construction, BASIL analysis, CBF, map smoothing, CVR registration, calculation and saving
clear all
clc

SUBJECT.masterdir='/Fridge/users/jeroen/MOYAMOYA/';
SUBJECT.locationBASILinfo=[SUBJECT.masterdir 'BASIL_OPTIONS.txt']; % location .txt file with addition model options for CBF quantification BASIL

SUBJECT.tau = 2; % Label duration
SUBJECT.N_BS = 4; % Number of background suppression pulses
SUBJECT.labeleff = 0.85; %PCASL label efficiency
SUBJECT.lambda = 0.9;%water partition fraction, Alsop MRM 2014
SUBJECT.T1t = 1.3;    %s T1 of tissue DEFAULT is 1.3 at 3T, Alsop MRM 2014
SUBJECT.T1b = 1.65;   %s T1 of arterial blood DEFAULT is 1.65 at 3T, Alsop MRM 2014
SUBJECT.LookLocker_correction_factor_perPLD = [0.4226, 0.3830, 0.3471, 0.3146, 0.2851]; % computed in LookLocker_PCASL_Simulation_Extended.m, in/matlab/Simulations/ASL folder, deltaMxy_ratio_LL_noLL_forBASIL
SUBJECT.FWHM = 6; % smoothing kernel size 6mm FWHM, for CBF, AAT
SUBJECT.FWHM_CVR = 8; % smoothing kernel size 8mm FWHM, for CVR maps
SUBJECT.FWHM_M0 = 4.5; % smoothing kernel size 4.5mm FWHM, for M0_forQCBF for manual quantification
SUBJECT.outlierFactor = 2.5; % outlierFactor x stds from mean CBF are considered outliers in step1 of Dolui et al SCORE method
SUBJECT.range_adult_cbf = [0 75]; % colourbar range for adult CBF values
SUBJECT.range_child_cbf = [0 125]; % colourbar range for child CBF values
SUBJECT.range_cvr = [-50 50]; % colourbar range for CVR values
SUBJECT.range_AAT = [0.5 2.5]; % time (s), arterial arrival time
SUBJECT.range_AATdelta = [-0.7 0.7];% delta time (s), delta arterial arrival time, postACZ - preACZ
SUBJECT.range_aCBV = [0 2]; % arterial blodo volume estimate in volume fraction voxel (%)

%%%%%%%%%%%%%%%%%%%%%%%%% 1. Subject information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get subject folder name, select folder containing all patient data

SUBJECT.SUBJECTdir = uigetdir(SUBJECT.masterdir,'Select subject folder');

% create folder paths
SUBJECT.ANATOMYdir = [SUBJECT.SUBJECTdir,'/ANATOMY/']; % T1 anatomy path
SUBJECT.MNIdir = [SUBJECT.masterdir 'MNI/']; % MNI path
SUBJECT.SUBJECTMNIdir = [SUBJECT.SUBJECTdir '/MNI/']; % MNI path
SUBJECT.DICOMdir = [SUBJECT.SUBJECTdir,'/DICOM/']; % DICOM  path
SUBJECT.NIFTIdir = [SUBJECT.SUBJECTdir,'/NIFTI/']; % NIFTI  path
SUBJECT.ASLdir = [SUBJECT.SUBJECTdir,'/ASL/']; % ASL path
SUBJECT.RESULTSdir = [SUBJECT.SUBJECTdir,'/RESULTS/']; % RESULTS path

% create folders
if ~isfolder({SUBJECT.ANATOMYdir; SUBJECT.DICOMdir; SUBJECT.NIFTIdir; SUBJECT.ASLdir; SUBJECT.RESULTSdir})
    mkdir(SUBJECT.ANATOMYdir); % create Anatomy folder
    mkdir(SUBJECT.SUBJECTMNIdir); % create subject MNI folder
    mkdir(SUBJECT.DICOMdir); % create DICOM folder
    mkdir(SUBJECT.NIFTIdir); % create NIFTI folder
    mkdir(SUBJECT.ASLdir); % create ASL folder
    mkdir(SUBJECT.RESULTSdir); % create RESULTS folder
end

% rename DICOM files in DICOM folder, and convert to NIFTI
ASLConvertDICOMtoNIFTI(SUBJECT.DICOMdir, SUBJECT.NIFTIdir)

% Get ASL nifti filenames
% preACZ path
filepreACZ = dir([SUBJECT.NIFTIdir, '*preACZ*2.nii.gz']);
% postACZ path
filepostACZ = dir([SUBJECT.NIFTIdir, '*postACZ*2.nii.gz']);

if size(filepreACZ,1) || size(filepostACZ,1) > 1
    warning(' More than 1 preAZC or postACZ ASL dataset found -  !! taking the last scanned dataset... !!')
end

SUBJECT.preACZfilenameNIFTI = filepreACZ(end,1).name;
SUBJECT.preACZfilenameDCM = filepreACZ(end,1).name(1:end-7);
SUBJECT.postACZfilenameNIFTI = filepostACZ(end,1).name;
SUBJECT.postACZfilenameDCM = filepostACZ(end,1).name(1:end-7);

%%%%%%%%%%%%%%%%%%%%%%%%% 2. Extract DICOM information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch scan parameters
SUBJECT = ASLExtractParamsDICOM(SUBJECT, SUBJECT.preACZfilenameDCM);
SUBJECT.dummyfilenameSaveNII = [SUBJECT.NIFTIdir SUBJECT.preACZfilenameNIFTI]; % lication .nii.gz NIFTI to be used as dummy template for saving NII's in the tool

%%%%%%%%%%%%%%%%%%%%%%%%% 3. Modify NIFTI to correct names, correct Mz loss (small fip angle) using Look Locker correction %%%%%%%%%%%%%%%%%%%%%
% save per PLD and control and label volumes (interleaved), and save all ASL and M0 in struct SUBJECT
SUBJECT = ASLPrepareM0ASLData(SUBJECT, SUBJECT.preACZfilenameNIFTI, 'preACZ'); % preACZ
SUBJECT = ASLPrepareM0ASLData(SUBJECT, SUBJECT.postACZfilenameNIFTI, 'postACZ'); % postACZ

disp('DICOMs converted to NIFTI');

%%%%%%%%%%%%%%%%%%%%%%%%%% 4. Generate T1 from M0 , T1 Tissue segmentation and registration to T1 anatomy and MNI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Locate T1 anatomy NIFTI
SUBJECT.T1ANATfilenameNIFTI = dir([SUBJECT.NIFTIdir, '*T1*.nii*']); % find T1 anatomy NIFTI filename
SUBJECT.T1ANATfilenameNIFTI = SUBJECT.T1ANATfilenameNIFTI.name;

% T1w ANATOMY scan brain extraction and FSL FAST segmentation in GM, WM and CSF of  into ANATOMY dir
T1Processing(SUBJECT, SUBJECT.T1ANATfilenameNIFTI);

% create T1fromM0 and save M0 from ASL multiPLD data and tissue segmentation using FSL FAST
SUBJECT = ASLT1fromM0Processing(SUBJECT, 'preACZ');
SUBJECT = ASLT1fromM0Processing(SUBJECT, 'postACZ');

% Register T1 and tissue segmentations to ASL space
ASLT1Registration(SUBJECT,'preACZ');
ASLT1Registration(SUBJECT,'postACZ');

% Register MNI segmentations to ASL space using the T1fromM0
ASLMNIRegistration(SUBJECT,'preACZ');
ASLMNIRegistration(SUBJECT,'postACZ');
disp('T1 ASL MNI tissue segmentation and registration finished');

%%%%%%%%%%%%%%%%%%%%%%%%%% 5. Outlier identification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dolui et al. SCORE outlier method
SUBJECT = ASLOutlierRemoval(SUBJECT, 'preACZ');
SUBJECT = ASLOutlierRemoval(SUBJECT, 'postACZ');


%%%%%%%%%%%%%%%%%%%%%%%%%% 6. BASIL CBF Analysis for both Original and Outlier removed ASL data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Perform BASIL analysis for both original and outlier removed ASL data')

session = {'preACZ', 'postACZ'};
for i=1:length(session)
    prefix = char(session(i));
    % %%%% % all PLD for AAT (arterial arrival time map) % %%%%%
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_allPLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_allPLD'], [1:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_allPLD_label1label2_OR.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_allPLD_OR'], [1:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
    % %%%% % 2tolast PLD for CBF map % %%%% %
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_2tolastPLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_2tolastPLD_forCBF'], [2:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_2tolastPLD_label1label2_OR.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_2tolastPLD_forCBF_OR'], [2:SUBJECT.NPLDS], SUBJECT.locationBASILinfo)
    % %%%% % 1to2 PLDs for ATA map -> compute COV --> then do no fit for the arterial component 'artoff' % %%%% %    
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_1to2PLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_1to2PLD_artoff_forATA'], [1:2], SUBJECT.locationBASILinfo,'artoff') 
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_1to2PLD_label1label2_OR.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_1to2PLD_artoff_forATA_OR'], [1:2], SUBJECT.locationBASILinfo,'artoff')
    % %%%% % 1to2 PLDs for aCBV map  then do fit for the arterial component 'noartoff'%  %%%% %
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_1to2PLD_label1label2.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_1to2PLD_noartoff_foraCBV'], [1:2], SUBJECT.locationBASILinfo, 'noartoff')    
    ASLBASILanalysis(SUBJECT, [SUBJECT.ASLdir prefix '_1to2PLD_label1label2_OR.nii.gz'], [SUBJECT.ASLdir prefix '_M0.nii.gz'], [SUBJECT.ASLdir prefix '_BASIL_1to2PLD_noartoff_foraCBV_OR'], [1:2], SUBJECT.locationBASILinfo, 'noartoff')   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  7. Generate resulting CBF/CVR/AAT/aCBV .png images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJECT = ASLSaveResultsCBFAATCVR(SUBJECT, '');
SUBJECT = ASLSaveResultsCBFAATCVR(SUBJECT, '_OR'); % save outlier removed data

%%%%%%%%%%%%%%%%%%%%%%%%%% 8. Save workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Workspace_ClinicalASL = [SUBJECT.ASLdir '/Workspace_ClinicalASL.mat'];
save(Workspace_ClinicalASL,'-v7.3');
disp('-- Finished --');
toc
