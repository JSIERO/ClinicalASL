function SUBJECT = ASLT1fromM0Processing(SUBJECT, prefix)
% ClinicalASL toolbox 2023, JCWSiero
% Construct T1w image from multi PLD M0 image

% brain extraction on M0 image
eval(['!bet ' SUBJECT.ASLdir prefix '_M0.nii.gz  ' SUBJECT.ASLdir prefix '_M0_brain.nii.gz' ' -m -f 0.4 -g 0']) % extract brain mask of M0's
eval(['!fslmaths ' SUBJECT.ASLdir prefix '_M0_brain_mask -dilM ' SUBJECT.ASLdir prefix '_temp_M0_brain_mask_dilM']);
eval(['!fslmaths ' SUBJECT.ASLdir prefix '_M0_brain_mask -kernel 2D -ero ' SUBJECT.ASLdir prefix '_temp_M0_brain_mask_ero']);
eval(['!fslmaths ' SUBJECT.ASLdir prefix '_M0 -mul ' SUBJECT.ASLdir prefix '_temp_M0_brain_mask_ero -kernel 2D -dilD ' SUBJECT.ASLdir prefix '_temp_M0_dilD']);

NII = load_untouch_nii([SUBJECT.ASLdir prefix '_M0_brain_mask.nii.gz']);
SUBJECT.(prefix).brainmask = logical(NII.img);

SUBJECT.(prefix).nanmask = double(SUBJECT.preACZ.brainmask);
SUBJECT.(prefix).nanmask(SUBJECT.(prefix).nanmask ==0) = NaN;

NII = load_untouch_nii([SUBJECT.ASLdir prefix '_temp_M0_dilD.nii.gz']);
M0_dilD = double(NII.img);

% smooth M0 with SUBJECT.FWHM_M0 kernel using edge preserving smoothin by Susan (FSL)
SUBJECT.(prefix).M0_forQCBF = ASLSmoothImage(M0_dilD, 3, SUBJECT.FWHM_M0, SUBJECT.VOXELSIZE);

SaveDataNII(SUBJECT.(prefix).M0_forQCBF , [SUBJECT.ASLdir prefix '_M0_forQCBF.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [0 500],SUBJECT.TR) % save T1fromM0
eval(['!rm ' SUBJECT.ASLdir '*temp*']); % remove temp files

% compute T1w image from multi-delay PCASL (Look-Locker) M0 data across the different PLDs
for i = 1:SUBJECT.NPLDS
M0_allPLD_noLLcorr(:,:,:,i) = SUBJECT.(prefix).M0_allPLD(:,:,:,i) * SUBJECT.LookLocker_correction_factor_perPLD(i); % remove the Look Locker correction to compute the T1w profile
end

T1fromM0 = ASLT1fromM0Compute(M0_allPLD_noLLcorr, SUBJECT.(prefix).brainmask, SUBJECT.PLDS);

SaveDataNII(T1fromM0, [SUBJECT.ASLdir prefix '_T1fromM0.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [0 500],SUBJECT.TR) % save T1fromM0

% tissue segment T1fromM0 using FSL FAST and load into SUBEJCT struct
eval(['!fast -b -g -B ' [SUBJECT.ASLdir prefix '_T1fromM0.nii.gz'] ' ' [SUBJECT.ASLdir prefix '_T1fromM0.nii.gz']]) 
eval(['!fslmaths ' SUBJECT.ASLdir prefix '_T1fromM0_restore.nii.gz ' SUBJECT.ASLdir prefix '_T1fromM0.nii.gz']); % use the _restore data as the T1fromM0

NII = load_untouch_nii([SUBJECT.ASLdir prefix '_T1fromM0.nii.gz']);
SUBJECT.(prefix).T1fromM0 = double(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir prefix '_T1fromM0_seg_0.nii.gz']);
SUBJECT.(prefix).CSFmask = logical(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir prefix '_T1fromM0_seg_1.nii.gz']);
SUBJECT.(prefix).GMmask = logical(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir prefix '_T1fromM0_seg_2.nii.gz']);
SUBJECT.(prefix).WMmask = logical(NII.img);
end
