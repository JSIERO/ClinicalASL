function [SUBJECT] = ASLExtractParamsDICOM_vTR(SUBJECT, filename)
% ClinicalASL toolbox 2023, JCWSiero
tic
disp('  Reading DICOM info... ')
info = dicominfo([SUBJECT.DICOMdir filename]);
elapsedTime = toc;
disp(['..this took: '  num2str(round(elapsedTime)) ' s'])

%%%%% 3. Scan information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJECT.TR=info.PerFrameFunctionalGroupsSequence.Item_1.CardiacSynchronizationSequence.Item_1.RRIntervalTimeNominal/1e3;% TR in s
SUBJECT.TE=info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.EchoTime;% TE in ms
SUBJECT.NDYNS = info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.NumberOfTemporalPositions; % number of repeated volumes including M0
SUBJECT.NPLDS = SUBJECT.NDYNS;
SUBJECT.VOXELSIZE = info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.PixelSpacing; % voxel size in mm
SUBJECT.VOXELSIZE(3) = info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.SliceThickness; % slice thickness in mm
SUBJECT.NSLICES = info.Private_2001_1018; % number of slices
SUBJECT.FLIPANGLE = info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.FlipAngle; % flip angle in degrees
SUBJECT.alpha_inv = SUBJECT.labeleff; %label efficiency
SUBJECT.alpha_BS = 0.95^SUBJECT.N_BS; %BS inversion efficiency, taken from Philips recon control parameters% Garcia MRM 2005
SUBJECT.alpha = SUBJECT.alpha_inv*SUBJECT.alpha_BS;
SUBJECT.([filename '_DCMinfo']) = info; % store DICOM info in SUBJECT