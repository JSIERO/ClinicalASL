function [SUBJECT] = ASLExtractParamsPARREC(SUBJECT)
% ClinicalASL toolbox 2023, JCWSiero

%%%%% 3. Scan information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJECT.TR = 5.5;% TR in s
SUBJECT.TE = 10.8;% TE in ms
SUBJECT.NPLDS = 5; % Amount of PLD's
SUBJECT.NDYNS = 24; % number of repeated volumes including M0
SUBJECT.NREPEATS = SUBJECT.NDYNS - 1; % number of repeated volumes without M0
SUBJECT.VOXELSIZE = [3 3]; % voxel size in mm
SUBJECT.VOXELSIZE(3) = 7; % slice thickness in mm
SUBJECT.NSLICES = 16; % number of slices
SUBJECT.FLIPANGLE = 25; % flip angle in degrees

% fetch timing  for every PLD, and for every slice and label control, and volume
%for t=1:SUBJECT.NSLICES*SUBJECT.NPLDS*SUBJECT.NDYNS*2
%SUBJECT.frametimes(t)=info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(t)).Private_2005_140f.Item_1.TriggerTime;
%end

SUBJECT.PLDS = [1.05, 1.618, 2.187,	2.755, 3.324]; %  PLDs in s, frametimes array should be sorted by PLD
SUBJECT.TIS = SUBJECT.PLDS + SUBJECT.tau; %  TIs in s

%UniqueSortFrametimes = unique(sort(SUBJECT.frametimes)); % take sort and unique operator for order across slices
SUBJECT.slicetime = 31.2; % slicetime in ms

SUBJECT.alpha_inv = SUBJECT.labeleff; %label efficiency
SUBJECT.alpha_BS = 0.95^SUBJECT.N_BS; %BS inversion efficiency, taken from Philips recon control parameters% Garcia MRM 2005
SUBJECT.alpha = SUBJECT.alpha_inv*SUBJECT.alpha_BS;
SUBJECT.TR_M0=SUBJECT.tau + SUBJECT.PLDS; %TR of the M0
