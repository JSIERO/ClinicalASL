function SUBJECT = ASLSaveCBFCVRResults_prepostACZ(SUBJECT, ORprefix)
%% CBF using 2ndtolast PLD, compute CVR
% register postACZ to preACZ, save BASIl results: NIFTI and PNG - loops over smoothed BASIl data

% Registration postACZ to preACZ space using T1fromM0 AFNI 3dAllineate, 12 dof wsinc interpolation, brain masks as weight
preACZ_T1fromM0_path = [SUBJECT.ASLdir 'preACZ_T1fromM0.nii.gz'];
postACZ_T1fromM0_path = [SUBJECT.ASLdir 'postACZ_T1fromM0.nii.gz'];
postACZ_T1fromM0_2preACZ_path = [SUBJECT.ASLdir 'postACZ_T1fromM0_2preACZ.nii.gz'];
preACZ_mask_path = [SUBJECT.ASLdir 'preACZ_M0_brain_mask.nii.gz'];
postACZ_mask_path = [SUBJECT.ASLdir 'postACZ_M0_brain_mask.nii.gz'];
postACZ_T1fromM0_2preACZ_mat = [SUBJECT.ASLdir 'postACZ_T1fromM0_2preACZ.aff12.1D'];
% Registration postACZ to preACZ
disp('Registration postACZ to preACZ CBF data')
eval(['!3dAllineate -input ' postACZ_T1fromM0_path ' -base ' preACZ_T1fromM0_path ' -prefix ' postACZ_T1fromM0_2preACZ_path ' -cost lpa -interp cubic -final wsinc5 -onepass -warp affine_general -weight ' preACZ_mask_path ' -source_mask ' postACZ_mask_path ' -1Dmatrix_save ' postACZ_T1fromM0_2preACZ_mat ' -overwrite']);
eval(['!fslcpgeom ' preACZ_T1fromM0_path ' ' postACZ_T1fromM0_2preACZ_path ' -d']);

preACZ_CBF_path = [SUBJECT.ASLdir 'preACZ_BASIL_2tolastPLD_forCBF' ORprefix '/native_space/perfusion_calib.nii.gz'];
postACZ_CBF_path = [SUBJECT.ASLdir 'postACZ_BASIL_2tolastPLD_forCBF' ORprefix '/native_space/perfusion_calib.nii.gz'];
postACZ_CBF_2preACZ_path = [SUBJECT.ASLdir 'postACZ_CBF' ORprefix '_2preACZ.nii.gz'];
eval(['!3dAllineate -input ' postACZ_CBF_path ' -master ' preACZ_CBF_path ' -prefix ' postACZ_CBF_2preACZ_path ' -1Dmatrix_apply ' postACZ_T1fromM0_2preACZ_mat ' -final wsinc5 -floatize -overwrite']);
eval(['!fslcpgeom ' preACZ_CBF_path ' ' postACZ_CBF_2preACZ_path ' -d']);
% register postACZ brain mask to preACZ
postACZ_2preACZ_mask_path = [ SUBJECT.ASLdir 'postACZ_M0_brain_mask_2preACZ.nii.gz'];
eval(['!3dAllineate -input ' postACZ_mask_path ' -master ' preACZ_mask_path ' -prefix ' postACZ_2preACZ_mask_path ' -1Dmatrix_apply ' postACZ_T1fromM0_2preACZ_mat ' -final NN -overwrite']);
eval(['!fslcpgeom ' preACZ_mask_path ' ' postACZ_2preACZ_mask_path ' -d']);
disp('Registration finished..')

% load BASIL CBF results
NII = load_untouch_nii(preACZ_CBF_path); SUBJECT.preACZ.(['CBF' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_CBF_path); SUBJECT.postACZ.(['CBF' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_CBF_2preACZ_path); SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ']) = double(NII.img);

% load registered psotaCZ to preACZ brainmask
NII = load_untouch_nii(postACZ_2preACZ_mask_path); SUBJECT.postACZ.brainmask_2preACZ = double(NII.img);

SUBJECT.preACZ.nanmask = double(SUBJECT.preACZ.brainmask);
SUBJECT.preACZ.nanmask(SUBJECT.preACZ.nanmask==0) = NaN;
SUBJECT.postACZ.nanmask = double(SUBJECT.postACZ.brainmask);
SUBJECT.postACZ.nanmask(SUBJECT.postACZ.nanmask==0) = NaN;
SUBJECT.postACZ.nanmask_2preACZ = SUBJECT.postACZ.brainmask_2preACZ;
SUBJECT.postACZ.nanmask_2preACZ(SUBJECT.postACZ.nanmask_2preACZ==0) = NaN;

% make combined mask from registered pre/post ACZ
SUBJECT.nanmask_reg = double(logical(SUBJECT.preACZ.brainmask + SUBJECT.postACZ.brainmask_2preACZ)); %combine preACZ and postACZ mask
SUBJECT.nanmask_reg(SUBJECT.nanmask_reg==0) = NaN;

% Compute CVR
SUBJECT.(['CVR' ORprefix]) = SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ']) - SUBJECT.preACZ.(['CBF' ORprefix]);

% CBF Smoothing
SUBJECT.preACZ.(['CBF' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.preACZ.(['CBF' ORprefix]).*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['CBF' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['CBF' ORprefix]).*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ']).*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

SUBJECT.(['CVR' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(['CVR' ORprefix]).*SUBJECT.nanmask_reg, 2, SUBJECT.FWHM_CVR, SUBJECT.VOXELSIZE);

% Compute CVR in percentage
SUBJECT.(['CVR' ORprefix '_percentage_smth']) = SUBJECT.(['CVR' ORprefix '_smth']) ./SUBJECT.preACZ.(['CBF' ORprefix '_smth']) .*SUBJECT.nanmask_reg * 100;
% 
%%% save final CBF, CVR NIFTI and .PNGs
smoothloop = {'', '_smth'};
for i=1:length(smoothloop)
  smthprefix = char(smoothloop(i));
  SaveDataNII(SUBJECT.preACZ.(['CBF' ORprefix smthprefix]), [SUBJECT.ASLdir 'preACZ_CBF' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['CBF' ORprefix smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF' ORprefix '_2preACZ' smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);

  SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF' ORprefix smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF' ORprefix smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' ORprefix smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' ORprefix smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' ORprefix '_2preACZ' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' ORprefix '_2preACZ' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
  if strcmp(smthprefix,'_smth')
    SaveDataNII(SUBJECT.(['CVR' ORprefix '_smth']), [SUBJECT.ASLdir 'CVR' ORprefix '_smth.nii.gz'], SUBJECT.dummyfilenameSaveNII,1,[], SUBJECT.TR);
    SaveDataNII(SUBJECT.(['CVR' ORprefix '_percentage_smth']),[SUBJECT.ASLdir 'CVR' ORprefix '_percentage_smth.nii.gz'], SUBJECT.dummyfilenameSaveNII,1,[], SUBJECT.TR);

    SaveFIGUREtoPNG(SUBJECT.(['CVR' ORprefix '_smth']), SUBJECT.nanmask_reg, SUBJECT.range_cvr, SUBJECT.RESULTSdir, ['CVR' ORprefix '_smth'],'CVR', 'vik');
    SaveFIGUREtoPNG(SUBJECT.(['CVR' ORprefix '_percentage_smth']), SUBJECT.nanmask_reg, [-100 100], SUBJECT.RESULTSdir, ['CVR' ORprefix '_PERCENTAGE_smth'],'%','vik');
  end
end
disp('CBF, CVR Results: NIFTI and .PNGs created')

%% Arrival transit time (AAT), aCBV, arterial transit artefact (ATA) for spatial COV computation, allPLD for another CBF map using all PLD
% register postACZ to preACZ, save BASIl results: NIFTI and PNG

% %%%% % all PLD for 'allCBF' map (arterial arrival time) map
preACZ_allPLD_CBF_path = [SUBJECT.ASLdir 'preACZ_BASIL_allPLD' ORprefix '/native_space/perfusion_calib.nii.gz'];
postACZ_allPLD_CBF_path = [SUBJECT.ASLdir 'postACZ_BASIL_allPLD' ORprefix '/native_space/perfusion_calib.nii.gz'];
postACZ_allPLD_CBF_2preACZ_path = [SUBJECT.ASLdir 'postACZ_BASIL_allPLD' ORprefix '/native_space/perfusion_calib_2preACZ.nii.gz'];
% %%%% % all PLD for AAT (arterial arrival time) map
preACZ_allPLD_AAT_path = [SUBJECT.ASLdir 'preACZ_BASIL_allPLD' ORprefix '/native_space/arrival.nii.gz'];
postACZ_allPLD_AAT_path = [SUBJECT.ASLdir 'postACZ_BASIL_allPLD' ORprefix '/native_space/arrival.nii.gz'];
postACZ_allPLD_AAT_2preACZ_path = [SUBJECT.ASLdir 'postACZ_BASIL_allPLD' ORprefix '/native_space/arrival_2preACZ.nii.gz'];
% %%%% % 1to2 PLD for ATA map to compute spatial COV
preACZ_1to2PLD_ATA_path = [SUBJECT.ASLdir 'preACZ_BASIL_1to2PLD_artoff_forATA' ORprefix '/native_space/perfusion_calib.nii.gz'];
postACZ_1to2PLD_ATA_path = [SUBJECT.ASLdir 'postACZ_BASIL_1to2PLD_artoff_forATA' ORprefix '/native_space/perfusion_calib.nii.gz'];
postACZ_1to2PLD_ATA_2preACZpath = [SUBJECT.ASLdir 'postACZ_BASIL_1to2PLD_artoff_forATA' ORprefix '/native_space/perfusion_calib_2preACZ.nii.gz'];
% %%%% % 1to2 PLD for aCBV map
preACZ_1to2PLD_aCBV_path = [SUBJECT.ASLdir 'preACZ_BASIL_1to2PLD_noartoff_foraCBV' ORprefix '/native_space/aCBV_calib.nii.gz'];
postACZ_1to2PLD_aCBV_path = [SUBJECT.ASLdir 'postACZ_BASIL_1to2PLD_noartoff_foraCBV' ORprefix '/native_space/aCBV_calib.nii.gz'];
postACZ_1to2PLD_aCBV_2preACZpath = [SUBJECT.ASLdir 'postACZ_BASIL_1to2PLD_noartoff_foraCBV' ORprefix '/native_space/aCBV_calib_2preACZ.nii.gz'];

% Registration postACZ to preACZ
disp('Registration postACZ to preACZ AAT, ATA, aCBV data')
eval(['!3dAllineate -input ' postACZ_allPLD_CBF_path ' -master ' preACZ_allPLD_CBF_path ' -prefix ' postACZ_allPLD_CBF_2preACZ_path ' -1Dmatrix_apply ' postACZ_T1fromM0_2preACZ_mat ' -final wsinc5 -floatize -overwrite']);
eval(['!fslcpgeom ' preACZ_allPLD_CBF_path ' ' postACZ_allPLD_CBF_2preACZ_path ' -d']);

eval(['!3dAllineate -input ' postACZ_allPLD_AAT_path ' -master ' preACZ_allPLD_AAT_path ' -prefix ' postACZ_allPLD_AAT_2preACZ_path ' -1Dmatrix_apply ' postACZ_T1fromM0_2preACZ_mat ' -final wsinc5 -floatize -overwrite']);
eval(['!fslcpgeom ' preACZ_allPLD_AAT_path ' ' postACZ_allPLD_AAT_2preACZ_path ' -d']);

eval(['!3dAllineate -input ' postACZ_1to2PLD_ATA_path ' -master ' preACZ_1to2PLD_ATA_path ' -prefix ' postACZ_1to2PLD_ATA_2preACZpath ' -1Dmatrix_apply ' postACZ_T1fromM0_2preACZ_mat ' -final wsinc5 -floatize -overwrite']);
eval(['!fslcpgeom ' preACZ_1to2PLD_ATA_path ' ' postACZ_1to2PLD_ATA_2preACZpath ' -d']);

eval(['!3dAllineate -input ' postACZ_1to2PLD_aCBV_path ' -master ' preACZ_1to2PLD_aCBV_path ' -prefix ' postACZ_1to2PLD_aCBV_2preACZpath ' -1Dmatrix_apply ' postACZ_T1fromM0_2preACZ_mat ' -final wsinc5 -floatize -overwrite']);
eval(['!fslcpgeom ' preACZ_1to2PLD_aCBV_path ' ' postACZ_1to2PLD_aCBV_2preACZpath ' -d']);
disp('Registration finished')

% load BASIL and registered postACZ to preACZ data
NII = load_untouch_nii(preACZ_allPLD_CBF_path); SUBJECT.preACZ.(['CBF_allPLD' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_allPLD_CBF_path); SUBJECT.postACZ.(['CBF_allPLD' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_allPLD_CBF_2preACZ_path); SUBJECT.postACZ.(['CBF_allPLD' ORprefix '_2preACZ']) = double(NII.img);

NII = load_untouch_nii(preACZ_allPLD_AAT_path); SUBJECT.preACZ.(['AAT' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_allPLD_AAT_path); SUBJECT.postACZ.(['AAT' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_allPLD_AAT_2preACZ_path); SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ']) = double(NII.img);

NII = load_untouch_nii(preACZ_1to2PLD_ATA_path); SUBJECT.preACZ.(['ATA' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_1to2PLD_ATA_path); SUBJECT.postACZ.(['ATA' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_1to2PLD_ATA_2preACZpath); SUBJECT.postACZ.(['ATA' ORprefix '_2preACZ']) = double(NII.img);

NII = load_untouch_nii(preACZ_1to2PLD_aCBV_path); SUBJECT.preACZ.(['aCBV' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_1to2PLD_aCBV_path); SUBJECT.postACZ.(['aCBV' ORprefix]) = double(NII.img);
NII = load_untouch_nii(postACZ_1to2PLD_aCBV_path); SUBJECT.postACZ.(['aCBV' ORprefix '_2preACZ']) = double(NII.img);

% Smooth data
SUBJECT.preACZ.(['CBF_allPLD' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.preACZ.(['CBF_allPLD' ORprefix]).*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['CBF_allPLD' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['CBF_allPLD' ORprefix]).*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['CBF_allPLD' ORprefix '_2preACZ_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['CBF_allPLD' ORprefix '_2preACZ']).*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

SUBJECT.preACZ.(['AAT' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.preACZ.(['AAT' ORprefix]).*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['AAT' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['AAT' ORprefix]).*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ']).*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

SUBJECT.preACZ.(['ATA' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.preACZ.(['ATA' ORprefix]).*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['ATA' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['ATA' ORprefix]).*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['ATA' ORprefix '_2preACZ_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['ATA' ORprefix '_2preACZ']).*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

SUBJECT.preACZ.(['aCBV' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.preACZ.(['aCBV' ORprefix]).*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['aCBV' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['aCBV' ORprefix]).*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['aCBV' ORprefix '_2preACZ_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['aCBV' ORprefix '_2preACZ']).*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

% compute deltaAAT map to reflect changes in AAT
SUBJECT.(['AATdelta' ORprefix '_smth'])  = ASLSmoothImage((SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ_smth']) - SUBJECT.preACZ.(['AAT' ORprefix '_smth'])).*SUBJECT.nanmask_reg, 2, SUBJECT.FWHM_CVR, SUBJECT.VOXELSIZE);

%%% save final CBF, CVR NIFTI and .PNGs
smoothloop = {'', '_smth'};
for i=1:length(smoothloop)
  smthprefix = char(smoothloop(i));  
  SaveDataNII(SUBJECT.preACZ.(['CBF_allPLD' ORprefix smthprefix]), [SUBJECT.ASLdir 'preACZ_CBF_allPLD' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['CBF_allPLD' ORprefix smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF_allPLD' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['CBF_allPLD' ORprefix '_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF_allPLD' ORprefix '_2preACZ' smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.preACZ.(['AAT' ORprefix smthprefix]), [SUBJECT.ASLdir 'preACZ_AAT' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['AAT' ORprefix smthprefix]), [SUBJECT.ASLdir 'postACZ_AAT' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_AAT' ORprefix '_2preACZ' smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.preACZ.(['ATA' ORprefix smthprefix]), [SUBJECT.ASLdir 'preACZ_ATA' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['ATA' ORprefix smthprefix]), [SUBJECT.ASLdir 'postACZ_ATA' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['ATA' ORprefix '_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_ATA' ORprefix '_2preACZ' smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.preACZ.(['aCBV' ORprefix smthprefix]), [SUBJECT.ASLdir 'preACZ_aCBV' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['aCBV' ORprefix smthprefix]), [SUBJECT.ASLdir 'postACZ_aCBV' ORprefix smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['aCBV' ORprefix '_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_aCBV' ORprefix '_2preACZ' smthprefix '.nii.gz'], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);

  SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF_allPLD' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF_allPLD' ORprefix smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF_allPLD' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF_allPLD' ORprefix '_2preACZ' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF_allPLD' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF_allPLD' ORprefix smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF_allPLD' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF_allPLD' ORprefix '_2preACZ' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.preACZ.(['AAT' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['preACZ_AAT' ORprefix smthprefix], 'time', 'devon');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['postACZ_AAT' ORprefix '_2preACZ' smthprefix], 'time', 'devon');
  SaveFIGUREtoPNG(SUBJECT.preACZ.(['ATA' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['preACZ_ATA' ORprefix smthprefix], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['ATA' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_ATA' ORprefix '_2preACZ' smthprefix], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.preACZ.(['aCBV' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['preACZ_aCBV' ORprefix smthprefix], '%', 'turku');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['aCBV' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_aCBV' ORprefix '_2preACZ' smthprefix], '%', 'turku');

  if strcmp(smthprefix,'_smth')
    SaveDataNII(SUBJECT.(['AATdelta' ORprefix '_smth']), [SUBJECT.ASLdir 'AATdelta' ORprefix '_smth.nii.gz'], SUBJECT.dummyfilenameSaveNII,1,[], SUBJECT.TR);
    
    SaveFIGUREtoPNG(SUBJECT.(['AATdelta' ORprefix '_smth']), SUBJECT.nanmask_reg, SUBJECT.range_AATdelta, SUBJECT.RESULTSdir, ['AATdelta' ORprefix '_smth'],'time_delta', 'vik');
  end
end

