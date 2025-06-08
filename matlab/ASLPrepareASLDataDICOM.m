function SUBJECT = ASLPrepareASLDataDICOM(SUBJECT, filename, prefix , fast)
% ClinicalASL toolbox 2023, JCWSiero
% Prepare (multidelay) ASL data: interleave control and label files per PLD, M0 and perform Look Locker Correction

DATA = double(niftiread([SUBJECT.NIFTIdir filename]));
DIMS = size(DATA);

SUBJECT.(prefix).M0ASL_allPLD = zeros(DIMS(1),DIMS(2), DIMS(3),SUBJECT.NDYNS, SUBJECT.NPLDS, 2);
SUBJECT.(prefix).ASL_label1label2_allPLD = zeros(DIMS(1),DIMS(2), DIMS(3),SUBJECT.NREPEATS*2, SUBJECT.NPLDS);

% Split (multidelay) ASL data into separate label and control files per PLD, correct Mz loss due to small flip angle using Look Locker correction factor
for i = 1:SUBJECT.NPLDS
    SUBJECT.(prefix).M0ASL_allPLD(:,:,:,1:SUBJECT.NDYNS,i,1) = DATA(:,:,:,i:2*SUBJECT.NPLDS:SUBJECT.NPLDS*SUBJECT.NDYNS*2)./ SUBJECT.LookLocker_correction_factor_perPLD(i);% for label images per PLD per volume
    SUBJECT.(prefix).M0ASL_allPLD(:,:,:,1:SUBJECT.NDYNS,i,2) = DATA(:,:,:,i+SUBJECT.NPLDS:2*SUBJECT.NPLDS:SUBJECT.NPLDS*SUBJECT.NDYNS*2)./ SUBJECT.LookLocker_correction_factor_perPLD(i);% for control images per PLD per volume

    SaveDataNII(squeeze(SUBJECT.(prefix).M0ASL_allPLD(:,:,:,:,i,1)), [SUBJECT.ASLdir prefix '_PLD0' num2str(i) '_label1'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(squeeze(SUBJECT.(prefix).M0ASL_allPLD(:,:,:,:,i,2)), [SUBJECT.ASLdir prefix '_PLD0' num2str(i) '_label2'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
end

% Interleave control and label ASL data per PLD, M0 and perform Look Locker Correction
for i = 1:SUBJECT.NPLDS
    SUBJECT.(prefix).ASL_label1label2_allPLD(:,:,:,:,i) = ASLInterleaveControlTag(squeeze(SUBJECT.(prefix).M0ASL_allPLD(:,:,:,2:end,i,1)), squeeze(SUBJECT.(prefix).M0ASL_allPLD(:,:,:,2:end,i,2)));
    SaveDataNII(squeeze(SUBJECT.(prefix).ASL_label1label2_allPLD(:,:,:,:,i)), [SUBJECT.ASLdir prefix '_PLD0' num2str(i) '_label1label2'], SUBJECT.dummyfilenameSaveNII, 1,[], SUBJECT.TR)  % save M0 image , only first PLDs
end

disp('Saving ASL data interleaved label control: all PLDs for AAT maps')
SaveDataNII(reshape(SUBJECT.(prefix).ASL_label1label2_allPLD, DIMS(1),DIMS(2), DIMS(3), SUBJECT.NPLDS*SUBJECT.NREPEATS*2), [SUBJECT.ASLdir prefix '_allPLD_label1label2'], SUBJECT.dummyfilenameSaveNII, 1,[], SUBJECT.TR)  % save interleaved control label and for all PLDs as 4th dimension

disp('Saving ASL data interleaved label control: 1-to-2 PLDs for ATA maps')
SaveDataNII(reshape(SUBJECT.(prefix).ASL_label1label2_allPLD(:,:,:,:,1:2), DIMS(1),DIMS(2), DIMS(3), 2*SUBJECT.NREPEATS*2), [SUBJECT.ASLdir prefix '_1to2PLD_label1label2'], SUBJECT.dummyfilenameSaveNII, 1,[], SUBJECT.TR)  % save interleaved control label and for all PLDs as 4th dimension

disp('Saving ASL data interleaved label control: 2-to-last PLDs for CBF maps, "free" of arterial transit artefacts')
SaveDataNII(reshape(SUBJECT.(prefix).ASL_label1label2_allPLD(:,:,:,:,2:end), DIMS(1),DIMS(2), DIMS(3), (SUBJECT.NPLDS-1)*SUBJECT.NREPEATS*2), [SUBJECT.ASLdir prefix '_2tolastPLD_label1label2'], SUBJECT.dummyfilenameSaveNII, 1,[], SUBJECT.TR)  % save interleaved control label and for all PLDs as 4th dimension

% Construct M0 calibration image
SUBJECT.(prefix).M0_allPLD = squeeze(mean(SUBJECT.(prefix).M0ASL_allPLD(:,:,:,1,:,:),6)); % average all PLD M0 control label images, and save in SUBJECT struct
SUBJECT.(prefix).M0 = squeeze(SUBJECT.(prefix).M0_allPLD(:,:,:,1)); % save 1st PLD M0 as true M0, and save in SUBJECT struct
if strcmp(fast, 'fast')
    % skip saving M0 image all PLDs
else
    SaveDataNII(SUBJECT.(prefix).M0_allPLD, [SUBJECT.ASLdir prefix '_M0_allPLD'], SUBJECT.dummyfilenameSaveNII, 1,[], SUBJECT.TR)  % save M0 image all PLDs
end
disp('Saving M0 image')
SaveDataNII(SUBJECT.(prefix).M0, [SUBJECT.ASLdir prefix '_M0'], SUBJECT.dummyfilenameSaveNII, 1,[], SUBJECT.TR)  % save M0 image , only first PLDs
