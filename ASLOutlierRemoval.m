function SUBJECT = ASLOutlierRemoval(SUBJECT, prefix)
% ClinicalASL toolbox 2023, JCWSiero
tic
outputmap = [SUBJECT.ASLdir prefix '_BASIL_perDYNAMIC/'];
if ~exist(outputmap,'dir')
    mkdir(outputmap)
end

DIMS = size(SUBJECT.(prefix).ASL_label1label2_allPLD);

if ~exist([outputmap 'DYNAMIC_1'],'dir')
    % Loop: compute CBF per dynamics, making sure we have 2 (control,label) x nPLD sized datasets
    for i=1:SUBJECT.NREPEATS
        %save data for BASIL analysis
        SaveDataNII(reshape(SUBJECT.(prefix).ASL_label1label2_allPLD(:,:,:,2*i-1:2*i,:),DIMS(1),DIMS(2),DIMS(3), SUBJECT.NPLDS*2), [outputmap prefix '_oneDYNAMIC_label1label2.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [],SUBJECT.TR)

        % perform BASIL analysis
        locationASL_multiPLD = [outputmap prefix '_oneDYNAMIC_label1label2.nii.gz'];
        locationM0 = [SUBJECT.ASLdir prefix '_M0.nii.gz'];
        locationMask = [SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz'];
        disp(['****************************************************************************************  BASIL analysis on DYNAMIC = ' num2str(i) '  *****************************************************'])

        ASLBASILanalysis(SUBJECT, locationASL_multiPLD, locationM0, locationMask, [outputmap 'DYNAMIC_' num2str(i)], [1:SUBJECT.NPLDS], SUBJECT.locationBASILinfo, 'artoff')

        NII = load_untouch_nii([SUBJECT.ASLdir prefix '_BASIL_perDYNAMIC/DYNAMIC_' num2str(i) '/native_space/perfusion_calib.nii.gz']);
        SUBJECT.(prefix).CBF_DYNAMIC(:,:,:,i) = double(NII.img);
    end

    % save CBF  per dynamic
    SaveDataNII(SUBJECT.(prefix).CBF_DYNAMIC,[SUBJECT.ASLdir prefix '_CBF_perDYNAMIC.nii.gz'], SUBJECT.dummyfilenameSaveNII ,1,[],SUBJECT.TR);
    disp(['Finished saving CBF data per dynamic: ' prefix]);
else
    disp(['Seems BASIL analysis per Dynamic has already been done, going straight to performing ASL oulierremoval by Duloi et al...' newline])
    disp(['Loading previous BASIL analysis per dynamic...'])
    for i=1:SUBJECT.NREPEATS
        NII = load_untouch_nii([SUBJECT.ASLdir prefix '_BASIL_perDYNAMIC/DYNAMIC_' num2str(i) '/native_space/perfusion_calib.nii.gz']);
        SUBJECT.(prefix).CBF_DYNAMIC(:,:,:,i) = double(NII.img);
    end
end

% %% ASL OutlierRemoval (volumes), SCORE method by Duloi et al JMRI 2017 %%%
[IOR_allsteps, IOR_step1, IOR_step2, NoOR_logical] = ASLOutlierRemovalPerform(SUBJECT.(prefix).CBF_DYNAMIC, SUBJECT.(prefix).brainmask, SUBJECT.outlierFactor, SUBJECT.(prefix).GMmask, SUBJECT.(prefix).WMmask, SUBJECT.(prefix).CSFmask);

if isempty(IOR_allsteps)
    IOR_allsteps=0;
end
% Write identified outliers to .txt
writematrix([numel(IOR_step1), numel(IOR_step2)],[SUBJECT.ASLdir prefix '_OutlierInformation_nvols_step1step2.txt'],'delimiter',' ');
writematrix(IOR_allsteps,[SUBJECT.ASLdir prefix '_OutlierInformation_whichvols_step1step2.txt'],'delimiter',' ');

%%%%%%%%%%%%%%%%%%%%%%%%%% Outlier removal and save to SUBJECT struct %%%%%%%%%%%%%%%%%%%%%%
SUBJECT.(prefix).NoOutliers_logical = NoOR_logical;
SUBJECT.(prefix).Ioutlier_allsteps = IOR_allsteps;
SUBJECT.(prefix).NREPEATS_OR = sum(SUBJECT.(prefix).NoOutliers_logical);
SUBJECT.(prefix).CBF_DYNAMIC_OR = SUBJECT.(prefix).CBF_DYNAMIC(:,:,:,SUBJECT.(prefix).NoOutliers_logical);
% save CBF outlier removed per dynamic
SaveDataNII(SUBJECT.(prefix).CBF_DYNAMIC_OR,[SUBJECT.ASLdir prefix '_CBF_perDYNAMIC_OR.nii.gz'], SUBJECT.dummyfilenameSaveNII ,1,[],SUBJECT.TR);

% remove outliers
NoOutliers_logical_forinterleaved = logical(zeros(1,2*SUBJECT.NREPEATS));
NoOutliers_logical_forinterleaved(1:2:end)=SUBJECT.(prefix).NoOutliers_logical; % odd indices
NoOutliers_logical_forinterleaved(2:2:end)=SUBJECT.(prefix).NoOutliers_logical; % even indices

SUBJECT.(prefix).ASL_label1label2_allPLD_OR = SUBJECT.(prefix).ASL_label1label2_allPLD(:,:,:,NoOutliers_logical_forinterleaved,:);

disp('Saving ASL data interleaved label control - OUTLIER REMOVED: all PLDs')
SaveDataNII(reshape(SUBJECT.(prefix).ASL_label1label2_allPLD_OR, DIMS(1),DIMS(2), DIMS(3), SUBJECT.NPLDS*SUBJECT.(prefix).NREPEATS_OR*2), [SUBJECT.ASLdir prefix '_allPLD_label1label2_OR.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1,[], SUBJECT.TR)  % save interleaved control label and for all PLDs as 4th dimension

disp('Saving ASL data interleaved label control - OUTLIER REMOVED: 2-to-last PLDs for CBF maps, "free" of arterial transit artefacts')
SaveDataNII(reshape(SUBJECT.(prefix).ASL_label1label2_allPLD_OR(:,:,:,:,2:end), DIMS(1),DIMS(2), DIMS(3), (SUBJECT.NPLDS-1)*SUBJECT.(prefix).NREPEATS_OR*2), [SUBJECT.ASLdir prefix '_2tolastPLD_label1label2_OR.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1,[], SUBJECT.TR)  % save interleaved control label and for all PLDs as 4th dimension

disp('Saving ASL data interleaved label control - OUTLIER REMOVED: 1-to-2 PLDs  for ATA, aCBF maps, and spatial COV ')
SaveDataNII(reshape(SUBJECT.(prefix).ASL_label1label2_allPLD_OR(:,:,:,:,1:2), DIMS(1),DIMS(2), DIMS(3), 2*SUBJECT.(prefix).NREPEATS_OR*2), [SUBJECT.ASLdir prefix '_1to2PLD_label1label2_OR.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1,[], SUBJECT.TR)  % save interleaved control label and for all PLDs as 4th dimension
toc
end
