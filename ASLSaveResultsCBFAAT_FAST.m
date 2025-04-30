function SUBJECT = ASLSaveResultsCBFAAT_FAST(SUBJECT, prefix)
% ClinicalASL toolbox 2023, JCWSiero
% save BASIl results: NIFTI and PNG - loops over smoothed BASIl data
% load BASIL output data : %% CBF Arrival transit time (AAT), aCBV, arterial transit artefact (ATA) for spatial COV computation, allPLD for another CBF map using all PLD

SUBJECT.(prefix).CBF = double(niftiread([SUBJECT.ASLdir prefix '_BASIL_2tolastPLD_forCBF/native_space/perfusion_calib.nii.gz']));
SUBJECT.(prefix).AAT = double(niftiread([SUBJECT.ASLdir prefix '_BASIL_allPLD_forAAT/native_space/arrival.nii.gz']));

SUBJECT.(prefix).nanmask = double(SUBJECT.(prefix).brainmask);
SUBJECT.(prefix).nanmask(SUBJECT.(prefix).nanmask==0) = NaN;

% Smooth data
SUBJECT.(prefix).AAT_smth = ASLSmoothImage(SUBJECT.(prefix).AAT.*SUBJECT.(prefix).nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

%%% save final CBF, AAT NIFTI and .PNGs

SaveDataNII(SUBJECT.(prefix).CBF, [SUBJECT.ASLdir 'CBF'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.(prefix).AAT, [SUBJECT.ASLdir 'AAT'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.(prefix).AAT_smth, [SUBJECT.ASLdir 'AAT_smth'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);

SaveFIGUREtoPNG(SUBJECT.(prefix).CBF, SUBJECT.(prefix).nanmask, SUBJECT.range_cbf, SUBJECT.RESULTSdir, ['CBF_' num2str(SUBJECT.range_cbf(2))], 'CBF', 'viridis');
SaveFIGUREtoPNG(SUBJECT.(prefix).AAT, SUBJECT.(prefix).nanmask, SUBJECT.range_AAT, SUBJECT.RESULTSdir, 'AAT', 'time', 'devon');
SaveFIGUREtoPNG(SUBJECT.(prefix).AAT_smth, SUBJECT.(prefix).nanmask, SUBJECT.range_AAT, SUBJECT.RESULTSdir, 'AAT_smth', 'time', 'devon');

disp('CBF, AAT Results: NIFTI and .PNGs created')
