function SUBJECT = ASLSaveResultsCBFAAT(SUBJECT, prefix, ORprefix)
% ClinicalASL toolbox 2023, JCWSiero
% save BASIl results: NIFTI and PNG - loops over smoothed BASIl data
% load BASIL output data : %% CBF Arrival transit time (AAT), aCBV, arterial transit artefact (ATA) for spatial COV computation, allPLD for another CBF map using all PLD

NII = load_untouch_nii([SUBJECT.ASLdir prefix '_BASIL_2tolastPLD_forCBF' ORprefix '/native_space/perfusion_calib.nii.gz']); SUBJECT.(prefix).(['CBF' ORprefix]) = double(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir prefix '_BASIL_allPLD_forAAT' ORprefix '/native_space/perfusion_calib.nii.gz']); SUBJECT.(prefix).(['CBF_allPLD' ORprefix]) = double(NII.img); % also procude CBF map using all PLDs
NII = load_untouch_nii([SUBJECT.ASLdir prefix '_BASIL_allPLD_forAAT' ORprefix '/native_space/arrival.nii.gz']); SUBJECT.(prefix).(['AAT' ORprefix]) = double(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir prefix '_BASIL_1to2PLD_forATA' ORprefix '/native_space/perfusion_calib.nii.gz']); SUBJECT.(prefix).(['ATA' ORprefix]) = double(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir prefix '_BASIL_1to2PLD_foraCBV' ORprefix '/native_space/aCBV_calib.nii.gz']); SUBJECT.(prefix).(['aCBV' ORprefix]) = double(NII.img);

SUBJECT.(prefix).nanmask = double(SUBJECT.(prefix).brainmask);
SUBJECT.(prefix).nanmask(SUBJECT.(prefix).nanmask==0) = NaN;

% Smooth data
SUBJECT.(prefix).(['CBF' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(prefix).(['CBF' ORprefix]).*SUBJECT.(prefix).nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.(prefix).(['CBF_allPLD' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(prefix).(['CBF_allPLD' ORprefix]).*SUBJECT.(prefix).nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.(prefix).(['AAT' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(prefix).(['AAT' ORprefix]).*SUBJECT.(prefix).nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.(prefix).(['ATA' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(prefix).(['ATA' ORprefix]).*SUBJECT.(prefix).nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.(prefix).(['aCBV' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(prefix).(['aCBV' ORprefix]).*SUBJECT.(prefix).nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth


%%% save final CBF, CVR NIFTI and .PNGs
smoothloop = {'', '_smth'};
for  i=1:length(smoothloop)
    smthprefix = char(smoothloop(i));
    SaveDataNII(SUBJECT.(prefix).(['CBF' ORprefix smthprefix]), [SUBJECT.ASLdir 'CBF' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.(prefix).(['CBF_allPLD' ORprefix smthprefix]), [SUBJECT.ASLdir 'CBF_allPLD' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.(prefix).(['AAT' ORprefix smthprefix]), [SUBJECT.ASLdir 'AAT' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.(prefix).(['ATA' ORprefix smthprefix]), [SUBJECT.ASLdir 'ATA' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.(prefix).(['aCBV' ORprefix smthprefix]), [SUBJECT.ASLdir 'aCBV' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    
    SaveFIGUREtoPNG(SUBJECT.(prefix).(['CBF' ORprefix smthprefix]), SUBJECT.(prefix).nanmask, SUBJECT.range_cbf, SUBJECT.RESULTSdir, ['CBF' ORprefix smthprefix '_' num2str(SUBJECT.range_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.(prefix).(['CBF_allPLD' ORprefix smthprefix]), SUBJECT.(prefix).nanmask, SUBJECT.range_cbf, SUBJECT.RESULTSdir, ['CBF_allPLD' ORprefix smthprefix '_' num2str(SUBJECT.range_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.(prefix).(['AAT' ORprefix smthprefix]), SUBJECT.(prefix).nanmask, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['AAT'  ORprefix smthprefix], 'time', 'devon');
    SaveFIGUREtoPNG(SUBJECT.(prefix).(['ATA' ORprefix smthprefix]), SUBJECT.(prefix).nanmask, SUBJECT.range_cbf, SUBJECT.RESULTSdir, ['ATA' ORprefix smthprefix], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.(prefix).(['aCBV' ORprefix smthprefix]), SUBJECT.(prefix).nanmask, SUBJECT.range_aCBV, SUBJECT.RESULTSdir, ['aCBV' ORprefix smthprefix], '%', 'turku');

end
disp('CBF, AAT, ATA, aCBV Results: NIFTI and .PNGs created')
