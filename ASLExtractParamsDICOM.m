function [SUBJECT] = ASLExtractParamsDICOM(SUBJECT, filename)
% ClinicalASL toolbox 2023, JCWSiero
tic
disp('  Reading DICOM info... ')
info = dicominfo([SUBJECT.DICOMdir filename]);
elapsedTime = toc;
disp(['..this took: '  num2str(round(elapsedTime)) ' s'])

%%%%% 3. Scan information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJECT.TR=info.PerFrameFunctionalGroupsSequence.Item_1.CardiacSynchronizationSequence.Item_1.RRIntervalTimeNominal/1e3;% TR in s
SUBJECT.TE=info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.EchoTime;% TE in ms
SUBJECT.NPLDS = info.Private_2001_1017; % Amount of PLD's
SUBJECT.NDYNS = info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.NumberOfTemporalPositions; % number of repeated volumes including M0
SUBJECT.NREPEATS = SUBJECT.NDYNS - 1; % number of repeated volumes without M0
SUBJECT.VOXELSIZE = info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.PixelSpacing; % voxel size in mm
SUBJECT.VOXELSIZE(3) = info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.SliceThickness; % slice thickness in mm
SUBJECT.NSLICES = info.Private_2001_1018; % number of slices
SUBJECT.FLIPANGLE = info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.FlipAngle; % flip angle in degrees

% fetch timing  for every PLD, and for every slice and label control, and volume
for t=1:SUBJECT.NSLICES*SUBJECT.NPLDS*SUBJECT.NDYNS*2
SUBJECT.frametimes(t)=info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(t)).Private_2005_140f.Item_1.TriggerTime;
end
SUBJECT.PLDS = SUBJECT.frametimes(1:SUBJECT.NPLDS)/1e3; %  PLDs in s, frametimes array should be sorted by PLD
SUBJECT.TIS = SUBJECT.PLDS + SUBJECT.tau; %  TIs in s

UniqueSortFrametimes = unique(sort(SUBJECT.frametimes)); % take sort and unique operator for order across slices
SUBJECT.slicetime = mean(diff(UniqueSortFrametimes(1:SUBJECT.NSLICES))); % slicetime in ms

SUBJECT.alpha_inv = SUBJECT.labeleff; %label efficiency
SUBJECT.alpha_BS = 0.95^SUBJECT.N_BS; %BS inversion efficiency, taken from Philips recon control parameters% Garcia MRM 2005
SUBJECT.alpha = SUBJECT.alpha_inv*SUBJECT.alpha_BS;
SUBJECT.TR_M0=SUBJECT.tau + SUBJECT.PLDS; %TR of the M0
SUBJECT.([filename '_DCMinfo']) = info; % store DICOM info in SUBJdCT