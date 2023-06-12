function SUBJECT = ASLSaveCBFCVRResults(SUBJECT, ORprefix)
% save BASIl results: NIFTI and PNG - loops over smoothed BASIl data
% load BASIL output data : %% CBF Arrival transit time (AAT), aCBV, arterial transit artefact (ATA) for spatial COV computation, allPLD for another CBF map using all PLD

NII = load_untouch_nii([SUBJECT.ASLdir 'BASIL' ORprefix '/native_space/perfusion_calib.nii.gz']); SUBJECT.(['CBF' ORprefix]) = double(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir 'BASIL_allPLD' ORprefix '/native_space/perfusion_calib.nii.gz']); SUBJECT.(['CBF_allPLD' ORprefix]) = double(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir 'BASIL_allPLD' ORprefix '/native_space/arrival.nii.gz']); SUBJECT.(['AAT' ORprefix]) = double(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir 'BASIL_1to2PLD' ORprefix '/native_space/perfusion_calib.nii.gz']); SUBJECT.(['ATA' ORprefix]) = double(NII.img);
NII = load_untouch_nii([SUBJECT.ASLdir 'BASIL_1to2PLD' ORprefix '/native_space/aCBV_calib.nii.gz']); SUBJECT.(['aCBV' ORprefix]) = double(NII.img);

% Smooth data
SUBJECT.(['CBF' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(['CBF' ORprefix]).*SUBJECT.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.(['CBF_allPLD' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(['CBF_allPLD' ORprefix]).*SUBJECT.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.(['AAT' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(['AAT' ORprefix]).*SUBJECT.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.(['ATA' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(['ATA' ORprefix]).*SUBJECT.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.(['aCBV' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(['aCBV' ORprefix]).*SUBJECT.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

SUBJECT.nanmask = double(SUBJECT.brainmask);
SUBJECT.nanmask(SUBJECT.nanmask==0) = NaN;

%%% save final CBF, CVR NIFTI and .PNGs
smoothloop = {'', '_smth'};
for  i=1:length(smoothloop)
    smthprefix = char(smoothloop(i));
    SaveDataNII(SUBJECT.(['CBF' ORprefix smthprefix]), [SUBJECT.ASLdir 'CBF' smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.(['CBF_allPLD' ORprefix smthprefix]), [SUBJECT.ASLdir 'CBF_allPLD' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.(['AAT' ORprefix smthprefix]), [SUBJECT.ASLdir 'AAT' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.(['ATA' ORprefix smthprefix]), [SUBJECT.ASLdir 'ATA' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.(['aCBV' ORprefix smthprefix]), [SUBJECT.ASLdir 'aCBV' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    
    
    SaveFIGUREtoPNG(SUBJECT.(['CBF' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['CBF' ORprefix smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.(['CBF' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['CBF' ORprefix smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');    
    SaveFIGUREtoPNG(SUBJECT.(['CBF_allPLD' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['CBF_allPLD' ORprefix smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.(['CBF_allPLD' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['CBF_allPLD' ORprefix smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.(['AAT' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['AAT'  ORprefix smthprefix], 'time', 'devon');
    SaveFIGUREtoPNG(SUBJECT.(['ATA' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['ATA' ORprefix smthprefix], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.(['aCBV' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['aCBV' ORprefix smthprefix], '%', 'turku');

end
disp('CBF, AAT, ATA, aCBV Results: NIFTI and .PNGs created')
