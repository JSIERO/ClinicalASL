function SUBJECT = ASLSaveResultsCBFAATCVR_vTRASL(SUBJECT)
% ClinicalASL toolbox 2023, JCWSiero

% copy NIFTIs to ASL_vTR folder
% M0
system(['cp ' SUBJECT.NIFTIdir '/' SUBJECT.preACZM0filenameNIFTI ' ' SUBJECT.ASLdir '/preACZ_M0.nii.gz']);
system(['cp ' SUBJECT.NIFTIdir '/' SUBJECT.postACZM0filenameNIFTI ' ' SUBJECT.ASLdir '/postACZ_M0.nii.gz']);
system(['fslroi ' SUBJECT.ASLdir '/preACZ_M0.nii.gz' ' ' SUBJECT.ASLdir '/preACZ_M0.nii.gz 0 1']);
system(['fslroi ' SUBJECT.ASLdir '/postACZ_M0.nii.gz' ' ' SUBJECT.ASLdir '/postACZ_M0.nii.gz 0 1']);
preACZ_M0_path = [SUBJECT.ASLdir '/preACZ_M0.nii.gz'];
postACZ_M0_path = [SUBJECT.ASLdir '/postACZ_M0.nii.gz'];
% AAT
SUBJECT.preACZAATfilenameNIFTI = [SUBJECT.preACZfilenameNIFTI(1:end-7) 'a.nii.gz'];
SUBJECT.postACZAATfilenameNIFTI = [SUBJECT.postACZfilenameNIFTI(1:end-7) 'a.nii.gz'];
system(['fslmaths ' SUBJECT.NIFTIdir '/' SUBJECT.preACZAATfilenameNIFTI ' -div 1000 ' SUBJECT.ASLdir '/preACZ_AAT.nii.gz']); % convert to s
system(['fslmaths ' SUBJECT.NIFTIdir '/' SUBJECT.postACZAATfilenameNIFTI ' -div 1000 ' SUBJECT.ASLdir '/postACZ_AAT.nii.gz']);
% CBF
SUBJECT.preACZCBFfilenameNIFTI = [SUBJECT.preACZfilenameNIFTI(1:end-7) '.nii.gz'];
SUBJECT.postACZCBFfilenameNIFTI = [SUBJECT.postACZfilenameNIFTI(1:end-7) '.nii.gz'];
system(['fslroi ' SUBJECT.NIFTIdir '/' SUBJECT.preACZfilenameNIFTI ' ' SUBJECT.ASLdir '/preACZ_CBF.nii.gz' ' ' '1 1']);
system(['fslroi ' SUBJECT.NIFTIdir '/' SUBJECT.postACZfilenameNIFTI ' ' SUBJECT.ASLdir '/postACZ_CBF.nii.gz' ' ' '1 1']);
preACZ_CBF_path = [SUBJECT.ASLdir 'preACZ_CBF.nii.gz'];
postACZ_CBF_path = [SUBJECT.ASLdir 'postACZ_CBF.nii.gz'];

% make brain masks
system(['bet2 ' SUBJECT.ASLdir '/preACZ_M0 ' SUBJECT.ASLdir '/preACZ_M0_brain -m -f 0.5']);
system(['bet2 ' SUBJECT.ASLdir '/postACZ_M0 ' SUBJECT.ASLdir '/postACZ_M0_brain  -m -f 0.5']);

preACZ_mask_path = [SUBJECT.ASLdir 'preACZ_M0_brain_mask.nii.gz'];
postACZ_mask_path = [SUBJECT.ASLdir 'postACZ_M0_brain_mask.nii.gz'];

%% compute CVR
% Registration postACZ to preACZ space using M0 AFNI 3dAllineate, 6 dof wsinc interpolation, brain masks as weight
postACZ_M0_2preACZ_path = [SUBJECT.ASLdir 'postACZ_M0_2preACZ.nii.gz'];
postACZ_M0_2preACZ_mat = [SUBJECT.ASLdir 'postACZ_M0_2preACZ.6dof.1D'];
postACZ_CBF_2preACZ_path = [SUBJECT.ASLdir 'postACZ_CBF_2preACZ.nii.gz'];
postACZ_mask_2preACZ_path = [SUBJECT.ASLdir 'postACZ_M0_brain_mask_2preACZ.nii.gz'];

% Registration postACZ to preACZ M0
disp('Registration postACZ to preACZ CBF data')
system(['3dAllineate -input ' postACZ_M0_path ' -base ' preACZ_M0_path ' -prefix ' postACZ_M0_2preACZ_path ' -cost lpa -autoweight -source_mask ' postACZ_mask_path ' -interp cubic -final wsinc5 -onepass -warp shift_rotate -1Dmatrix_save ' postACZ_M0_2preACZ_mat ' -overwrite']);
system(['fslcpgeom ' preACZ_M0_path ' ' postACZ_M0_2preACZ_path ' -d']);
SlicerPNGs(preACZ_M0_path, postACZ_M0_2preACZ_path, 'postACZ_M0', 'preACZ_M0', SUBJECT.ASLdir)
% Registration postACZ to preACZ CBF
system(['3dAllineate -input ' postACZ_CBF_path ' -master ' preACZ_CBF_path ' -prefix ' postACZ_CBF_2preACZ_path ' -1Dmatrix_apply ' postACZ_M0_2preACZ_mat ' -final wsinc5 -floatize -overwrite']);
system(['fslcpgeom ' preACZ_CBF_path ' ' postACZ_CBF_2preACZ_path ' -d']);
% register postACZ brain mask to preACZ
system(['3dAllineate -input ' postACZ_mask_path ' -master ' preACZ_mask_path ' -prefix ' postACZ_mask_2preACZ_path ' -1Dmatrix_apply ' postACZ_M0_2preACZ_mat ' -final NN -overwrite']);
system(['fslcpgeom ' preACZ_mask_path ' ' postACZ_2preACZ_mask_path ' -d']);
disp('Registration finished..')

% load CBF results
SUBJECT.preACZ.CBF = double(niftiread(preACZ_CBF_path));
SUBJECT.postACZ.CBF = double(niftiread(postACZ_CBF_path));
SUBJECT.postACZ.CBF_2preACZ = double(niftiread(postACZ_CBF_2preACZ_path));

% load brainmasks
SUBJECT.preACZ.brainmask = double(niftiread(preACZ_mask_path));
SUBJECT.postACZ.brainmask = double(niftiread(postACZ_mask_path));
SUBJECT.postACZ.brainmask_2preACZ = double(niftiread(postACZ_2preACZ_mask_path));

SUBJECT.preACZ.nanmask = double(SUBJECT.preACZ.brainmask);
SUBJECT.preACZ.nanmask(SUBJECT.preACZ.nanmask==0) = NaN;
SUBJECT.postACZ.nanmask = double(SUBJECT.postACZ.brainmask);
SUBJECT.postACZ.nanmask(SUBJECT.postACZ.nanmask==0) = NaN;
SUBJECT.postACZ.nanmask_2preACZ = SUBJECT.postACZ.brainmask_2preACZ;
SUBJECT.postACZ.nanmask_2preACZ(SUBJECT.postACZ.nanmask_2preACZ==0) = NaN;

% make combined mask from registered pre/post ACZ
SUBJECT.nanmask_reg = SUBJECT.preACZ.nanmask .*SUBJECT.postACZ.nanmask_2preACZ; %combine preACZ and postACZ mask

% Compute CVR
SUBJECT.CVR = SUBJECT.postACZ.CBF_2preACZ - SUBJECT.preACZ.CBF;

% CBF Smoothing
SUBJECT.preACZ.CBF_smth = ASLSmoothImage(SUBJECT.preACZ.CBF.*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.CBF_smth = ASLSmoothImage(SUBJECT.postACZ.CBF.*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.CBF_2preACZ_smth = ASLSmoothImage(SUBJECT.postACZ.CBF_2preACZ.*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

SUBJECT.CVR_smth = ASLSmoothImage(SUBJECT.CVR.*SUBJECT.nanmask_reg, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE);

% Compute CVR in percentage
SUBJECT.CVR_percentage_smth = SUBJECT.CVR_smth ./SUBJECT.preACZ.CBF_smth .*SUBJECT.nanmask_reg * 100;
%
%%% save final CBF, CVR NIFTI and .PNGs
smoothloop = {'', '_smth'};
for i=1:length(smoothloop)
  smthprefix = char(smoothloop(i));
  SaveDataNII(SUBJECT.preACZ.(['CBF' smthprefix]), [SUBJECT.ASLdir 'preACZ_CBF' smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['CBF' smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF' smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['CBF_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF_2preACZ' smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);

  SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF_2preACZ' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF_2preACZ' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
  if strcmp(smthprefix,'_smth')
    SaveDataNII(SUBJECT.CVR_smth, [SUBJECT.ASLdir 'CVR_smth'], SUBJECT.dummyfilenameSaveNII,1,[], SUBJECT.TR);
    SaveDataNII(SUBJECT.CVR_percentage_smth,[SUBJECT.ASLdir 'CVR_percentage_smth'], SUBJECT.dummyfilenameSaveNII,1,[], SUBJECT.TR);

    SaveFIGUREtoPNG(SUBJECT.CVR_smth, SUBJECT.nanmask_reg, SUBJECT.range_cvr, SUBJECT.RESULTSdir, 'CVR_smth','CVR', 'vik');
    SaveFIGUREtoPNG(SUBJECT.CVR_percentage_smth, SUBJECT.nanmask_reg, [-100 100], SUBJECT.RESULTSdir, 'CVR_PERCENTAGE_smth','%','vik');
  end
end
disp('CBF, CVR Results: NIFTI and .PNGs created')

%% Arrival transit time (AAT), 
% %%%% % all PLD for AAT (arterial arrival time) map
preACZ_AAT_path = [SUBJECT.ASLdir 'preACZ_AAT.nii.gz'];
postACZ_AAT_path = [SUBJECT.ASLdir 'postACZ_AAT.nii.gz'];
postACZ_AAT_2preACZ_path = [SUBJECT.ASLdir 'postACZ_AAT_2preACZ.nii.gz'];

% Registration postACZ to preACZ
disp('Registration postACZ to preACZ AAT,data')
system(['3dAllineate -input ' postACZ_AAT_path ' -master ' preACZ_AAT_path ' -prefix ' postACZ_AAT_2preACZ_path ' -1Dmatrix_apply ' postACZ_M0_2preACZ_mat ' -final wsinc5 -floatize -overwrite']);
system(['fslcpgeom ' preACZ_AAT_path ' ' postACZ_AAT_2preACZ_path ' -d']);
disp('Registration finished')

SUBJECT.preACZ.AAT = double(niftiread(preACZ_AAT_path));
SUBJECT.postACZ.AAT = double(niftiread(postACZ_AAT_path));
SUBJECT.postACZ.AAT_2preACZ = double(niftiread(postACZ_AAT_2preACZ_path));

% Smooth data
SUBJECT.preACZ.(['AAT' '_smth']) = ASLSmoothImage(SUBJECT.preACZ.AAT.*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['AAT' '_smth']) = ASLSmoothImage(SUBJECT.postACZ.AAT.*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['AAT' '_2preACZ_smth']) = ASLSmoothImage(SUBJECT.postACZ.AAT_2preACZ.*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

% compute deltaAAT map to reflect changes in AAT
SUBJECT.(['AATdelta' '_smth']) = ASLSmoothImage((SUBJECT.postACZ.(['AAT' '_2preACZ_smth']) - SUBJECT.preACZ.(['AAT' '_smth'])).*SUBJECT.nanmask_reg, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE);

%%% save final CBF, CVR NIFTI and .PNGs
smoothloop = {'', '_smth'};
for i=1:length(smoothloop)
  smthprefix = char(smoothloop(i));
  SaveDataNII(SUBJECT.preACZ.(['AAT' smthprefix]), [SUBJECT.ASLdir 'preACZ_AAT' smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['AAT' smthprefix]), [SUBJECT.ASLdir 'postACZ_AAT' smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
  SaveDataNII(SUBJECT.postACZ.(['AAT' '_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_AAT' '_2preACZ' smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR); 
  SaveFIGUREtoPNG(SUBJECT.preACZ.(['AAT' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['preACZ_AAT' smthprefix], 'time', 'devon');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['AAT' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['postACZ_AAT' smthprefix], 'time', 'devon');
  SaveFIGUREtoPNG(SUBJECT.postACZ.(['AAT' '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['postACZ_AAT' '_2preACZ' smthprefix], 'time', 'devon');
 
  if strcmp(smthprefix,'_smth')
    SaveDataNII(SUBJECT.(['AATdelta' '_smth']), [SUBJECT.ASLdir 'AATdelta' '_smth'], SUBJECT.dummyfilenameSaveNII,1,[], SUBJECT.TR); 
    SaveFIGUREtoPNG(SUBJECT.(['AATdelta' '_smth']), SUBJECT.nanmask_reg, SUBJECT.range_AATdelta, SUBJECT.RESULTSdir, ['AATdelta' '_smth'],'time_delta', 'vik');
  end
end
end

