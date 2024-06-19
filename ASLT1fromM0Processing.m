function SUBJECT = ASLT1fromM0Processing(SUBJECT, prefix)
% ClinicalASL toolbox 2023, JCWSiero
% Construct T1w image from multi PLD M0 image

% brain extraction on M0 image
system(['bet ' SUBJECT.ASLdir prefix '_M0  ' SUBJECT.ASLdir prefix '_M0_brain' ' -m -f 0.4 -g 0']); % extract brain mask of M0's
system(['fslmaths ' SUBJECT.ASLdir prefix '_M0_brain_mask -dilM ' SUBJECT.ASLdir prefix '_temp_M0_brain_mask_dilM']);
system(['fslmaths ' SUBJECT.ASLdir prefix '_M0_brain_mask -kernel 2D -ero ' SUBJECT.ASLdir prefix '_temp_M0_brain_mask_ero']);
system(['fslmaths ' SUBJECT.ASLdir prefix '_M0 -mul ' SUBJECT.ASLdir prefix '_temp_M0_brain_mask_ero -kernel 2D -dilD ' SUBJECT.ASLdir prefix '_temp_M0_dilD']);

SUBJECT.(prefix).brainmask = logical(niftiread([SUBJECT.ASLdir prefix '_M0_brain_mask']));

SUBJECT.(prefix).nanmask = double(SUBJECT.(prefix).brainmask);
SUBJECT.(prefix).nanmask(SUBJECT.(prefix).nanmask ==0) = NaN;

M0_dilD = double(niftiread([SUBJECT.ASLdir prefix '_temp_M0_dilD']));

% smooth M0 with SUBJECT.FWHM_M0 kernel using edge preserving smoothin by Susan (FSL)
SUBJECT.(prefix).M0_forQCBF = ASLSmoothImage(M0_dilD, 3, SUBJECT.FWHM_M0, SUBJECT.VOXELSIZE);

SaveDataNII(SUBJECT.(prefix).M0_forQCBF , [SUBJECT.ASLdir prefix '_M0_forQCBF'], SUBJECT.dummyfilenameSaveNII, 1, [0 500],SUBJECT.TR); % save T1fromM0
system(['rm ' SUBJECT.ASLdir '*temp*']); % remove temp files

% compute T1w image from multi-delay PCASL (Look-Locker) M0 data across the different PLDs
for i = 1:SUBJECT.NPLDS
M0_allPLD_noLLcorr(:,:,:,i) = SUBJECT.(prefix).M0_allPLD(:,:,:,i) * SUBJECT.LookLocker_correction_factor_perPLD(i); % remove the Look Locker correction to compute the T1w profile
end

T1fromM0 = ASLT1fromM0Compute(M0_allPLD_noLLcorr, SUBJECT.(prefix).brainmask, SUBJECT.PLDS);

SaveDataNII(T1fromM0, [SUBJECT.ASLdir prefix '_T1fromM0'], SUBJECT.dummyfilenameSaveNII, 1, [0 500],SUBJECT.TR); % save T1fromM0

% tissue segment T1fromM0 using FSL FAST and load into SUBEJCT struct
system(['fast -b -g -B ' [SUBJECT.ASLdir prefix '_T1fromM0.'] ' ' [SUBJECT.ASLdir prefix '_T1fromM0']]);
system(['fslmaths ' SUBJECT.ASLdir prefix '_T1fromM0_restore ' SUBJECT.ASLdir prefix '_T1fromM0']); % use the _restore data as the T1fromM0

SUBJECT.(prefix).T1fromM0  = double(niftiread([SUBJECT.ASLdir prefix '_T1fromM0']));
SUBJECT.(prefix).CSFmask = double(niftiread([SUBJECT.ASLdir prefix '_T1fromM0_seg_0']));
SUBJECT.(prefix).GMmask  = double(niftiread([SUBJECT.ASLdir prefix '_T1fromM0_seg_1']));
SUBJECT.(prefix).WMmask = double(niftiread([SUBJECT.ASLdir prefix '_T1fromM0_seg_2']));
end
