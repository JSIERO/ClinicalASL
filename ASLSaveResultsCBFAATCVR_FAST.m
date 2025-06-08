function SUBJECT = ASLSaveResultsCBFAATCVR_FAST(SUBJECT)
% ClinicalASL toolbox 2025, JCWSiero

% Registration postACZ to preACZ
if strcmp(SUBJECT.RegistrationMethod,'elastix')
    % use Elastix: S. Klein, M. Staring, K. Murphy, M.A. Viergever, J.P.W. Pluim, "elastix: a toolbox for intensity based medical image registration," IEEE Transactions on Medical Imaging, vol. 29, no. 1, pp. 196 - 205, January 2010

    disp('Registration T1fromM0 postACZ to preACZ data *********************************************************************')
    system(['elastix -f ' SUBJECT.preACZ_T1fromM0_path ' -m ' SUBJECT.postACZ_T1fromM0_path ' -fMask ' SUBJECT.preACZ_mask_path ' -p ' SUBJECT.ElastixParameterFile ' -loglevel info -out ' SUBJECT.ASLdir]);
    system(['mv -f ' fullfile(SUBJECT.ASLdir, 'result.0.nii.gz') ' ' SUBJECT.postACZ_T1fromM0_2preACZ_path]);
    system(['fslcpgeom ' SUBJECT.preACZ_T1fromM0_path ' ' SUBJECT.postACZ_T1fromM0_2preACZ_path ' -d']);

    disp('Registration CBF postACZ to preACZ *********************************************************************')
    system(['transformix -in ' SUBJECT.postACZ_CBF_path ' -out ' SUBJECT.ASLdir ' -tp ' fullfile(SUBJECT.ASLdir, 'TransformParameters.0.txt')]);
    system(['mv -f ' fullfile(SUBJECT.ASLdir, 'result.nii.gz') ' ' SUBJECT.postACZ_CBF_2preACZ_path]);
    system(['fslcpgeom ' SUBJECT.preACZ_CBF_path ' ' SUBJECT.postACZ_CBF_2preACZ_path ' -d']);

    disp('Registration AAT postACZ to preACZ *********************************************************************')
    system(['transformix -in ' SUBJECT.postACZ_AAT_path ' -out ' SUBJECT.ASLdir ' -tp ' fullfile(SUBJECT.ASLdir, 'TransformParameters.0.txt')]);
    system(['mv -f ' fullfile(SUBJECT.ASLdir, 'result.nii.gz') ' ' SUBJECT.postACZ_AAT_2preACZ_path]);
    system(['fslcpgeom ' SUBJECT.preACZ_AAT_path ' ' SUBJECT.postACZ_AAT_2preACZ_path ' -d']);
  
    disp('Registration ATA postACZ to preACZ *********************************************************************')
    system(['transformix -in ' SUBJECT.postACZ_ATA_path ' -out ' SUBJECT.ASLdir ' -tp ' fullfile(SUBJECT.ASLdir, 'TransformParameters.0.txt')]);
    system(['mv -f ' fullfile(SUBJECT.ASLdir, 'result.nii.gz') ' ' SUBJECT.postACZ_ATA_2preACZ_path]);
    system(['fslcpgeom ' SUBJECT.preACZ_ATA_path ' ' SUBJECT.postACZ_ATA_2preACZ_path ' -d']);

    disp('Registration mask postACZ to preACZ *********************************************************************')
    % change ResampleInterpolator to nearestneighbour for binary mask transformation to preACZ
    oldTextLine = '(ResampleInterpolator "FinalBSplineInterpolator")';
    newTextLine = '(ResampleInterpolator "FinalNearestNeighborInterpolator")';
    ChangeElastixParameterFileEntry(fullfile(SUBJECT.ASLdir, 'TransformParameters.0.txt'), oldTextLine, newTextLine, fullfile(SUBJECT.ASLdir, 'TransformParameters.0.NN.txt'));
    system(['transformix -in ' SUBJECT.postACZ_mask_path ' -out ' SUBJECT.ASLdir ' -tp ' fullfile(SUBJECT.ASLdir, 'TransformParameters.0.NN.txt')]);
    system(['mv -f ' fullfile(SUBJECT.ASLdir, 'result.nii.gz') ' ' SUBJECT.postACZ_mask_2preACZ_path]);
    system(['fslcpgeom ' SUBJECT.preACZ_mask_path ' ' SUBJECT.postACZ_mask_2preACZ_path ' -d']);
end
disp('Registration finished..')

% create figure .png to assess registration result
SlicerPNGs(SUBJECT.preACZ_T1fromM0_path, SUBJECT.postACZ_T1fromM0_2preACZ_path, 'postACZ_T1fromM0', 'preACZ_T1fromM0', SUBJECT.RESULTSdir)

% load BASIL CBF, AAT, masks and regsitration NIFTIs
SUBJECT.preACZ.CBF = double(niftiread(SUBJECT.preACZ_CBF_path));
SUBJECT.postACZ.CBF = double(niftiread(SUBJECT.postACZ_CBF_path));
SUBJECT.postACZ.CBF_2preACZ = double(niftiread(SUBJECT.postACZ_CBF_2preACZ_path));
SUBJECT.preACZ.AAT = double(niftiread(SUBJECT.preACZ_AAT_path));
SUBJECT.postACZ.AAT = double(niftiread(SUBJECT.postACZ_AAT_path));
SUBJECT.postACZ.AAT_2preACZ = double(niftiread(SUBJECT.postACZ_AAT_2preACZ_path));
SUBJECT.preACZ.ATA = double(niftiread(SUBJECT.preACZ_ATA_path));
SUBJECT.postACZ.ATA = double(niftiread(SUBJECT.postACZ_ATA_path));
SUBJECT.postACZ.ATA_2preACZ = double(niftiread(SUBJECT.postACZ_ATA_2preACZ_path));
SUBJECT.postACZ.brainmask_2preACZ = double(niftiread(SUBJECT.postACZ_mask_2preACZ_path));

% make combined mask from registered pre/post ACZ
SUBJECT.preACZ.nanmask = double(SUBJECT.preACZ.brainmask);
SUBJECT.preACZ.nanmask(SUBJECT.preACZ.nanmask==0) = NaN;
SUBJECT.postACZ.nanmask = double(SUBJECT.postACZ.brainmask);
SUBJECT.postACZ.nanmask(SUBJECT.postACZ.nanmask==0) = NaN;
SUBJECT.postACZ.nanmask_2preACZ = SUBJECT.postACZ.brainmask_2preACZ;
SUBJECT.postACZ.nanmask_2preACZ(SUBJECT.postACZ.nanmask_2preACZ==0) = NaN;
SUBJECT.nanmask_reg = SUBJECT.preACZ.nanmask .*SUBJECT.postACZ.nanmask_2preACZ; %combine preACZ and postACZ mask

%% Compute CVR
SUBJECT.CVR = SUBJECT.postACZ.CBF_2preACZ - SUBJECT.preACZ.CBF;

% Smooth AAT, CVR data
SUBJECT.CVR_smth = ASLSmoothImage(SUBJECT.CVR.*SUBJECT.nanmask_reg, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE);% 2D Smooth
SUBJECT.preACZ.AAT_smth = ASLSmoothImage(SUBJECT.preACZ.AAT.*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.AAT_smth = ASLSmoothImage(SUBJECT.postACZ.AAT.*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.AAT_2preACZ_smth = ASLSmoothImage(SUBJECT.postACZ.AAT_2preACZ.*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

%% save final CBF, AAT, CVR NIFTI and .PNGs

SaveDataNII(SUBJECT.preACZ.CBF, fullfile(SUBJECT.ASLdir, 'preACZ_CBF'), SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.postACZ.CBF, fullfile(SUBJECT.ASLdir, 'postACZ_CBF'), SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.postACZ.CBF_2preACZ, fullfile(SUBJECT.ASLdir, 'postACZ_CBF_2preACZ'), SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.preACZ.AAT, fullfile(SUBJECT.ASLdir, 'preACZ_AAT'), SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.postACZ.AAT, fullfile(SUBJECT.ASLdir, 'postACZ_AAT'), SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.postACZ.AAT_2preACZ, fullfile(SUBJECT.ASLdir, 'postACZ_AAT_2preACZ'), SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.preACZ.ATA, fullfile(SUBJECT.ASLdir, 'preACZ_ATA'), SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.postACZ.ATA, fullfile(SUBJECT.ASLdir, 'postACZ_ATA'), SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.postACZ.ATA_2preACZ, fullfile(SUBJECT.ASLdir, 'postACZ_ATA_2preACZ'), SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
SaveDataNII(SUBJECT.CVR_smth, fullfile(SUBJECT.ASLdir, 'CVR_smth'), SUBJECT.dummyfilenameSaveNII,1,[], SUBJECT.TR);

SaveFIGUREtoPNG(SUBJECT.preACZ.CBF, SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
SaveFIGUREtoPNG(SUBJECT.preACZ.CBF, SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
SaveFIGUREtoPNG(SUBJECT.postACZ.CBF_2preACZ, SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF_2preACZ_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
SaveFIGUREtoPNG(SUBJECT.postACZ.CBF_2preACZ, SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF_2preACZ_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
SaveFIGUREtoPNG(SUBJECT.CVR_smth, SUBJECT.nanmask_reg, SUBJECT.range_cvr, SUBJECT.RESULTSdir, 'CVR_smth','CVR', 'vik');
SaveFIGUREtoPNG(SUBJECT.preACZ.AAT_smth, SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, 'preACZ_AAT_smth', 'time', 'devon');
SaveFIGUREtoPNG(SUBJECT.postACZ.AAT_2preACZ_smth, SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, 'postACZ_AAT_2preACZ', 'time', 'devon');
SaveFIGUREtoPNG(SUBJECT.preACZ.ATA, SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, 'preACZ_ATA', 'CBF', 'viridis');
SaveFIGUREtoPNG(SUBJECT.postACZ.ATA_2preACZ, SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, 'postACZ_ATA_2preACZ', 'CBF', 'viridis');
disp('CBF, AAT, CVR Results: NIFTI and .PNGs created')

%% Save to DICOMS: CBF, AAT, ATA, CVR

% preACZ CBF
SaveDataDICOM(SUBJECT.preACZ.CBF, fullfile(SUBJECT.DICOMdir, SUBJECT.preACZfilenameDCM_CBF), fullfile(SUBJECT.DICOMRESULTSdir, [SUBJECT.preACZfilenameDCM_CBF '.dcm']), 'WIP preACZ CBF MD-ASL', SUBJECT.range_adult_cbf, 'CBF')
% preACZ AAT
SaveDataDICOM(SUBJECT.preACZ.AAT_smth, fullfile(SUBJECT.DICOMdir, SUBJECT.preACZfilenameDCM_AAT), fullfile(SUBJECT.DICOMRESULTSdir, [SUBJECT.preACZfilenameDCM_AAT '.dcm']), 'WIP preACZ AAT MD-ASL', SUBJECT.range_AAT, 'AAT')
% CVR
SaveDataDICOM(SUBJECT.CVR_smth, fullfile(SUBJECT.DICOMdir, SUBJECT.preACZfilenameDCM_CVR), fullfile(SUBJECT.DICOMRESULTSdir, [SUBJECT.preACZfilenameDCM_CVR '.dcm']), 'WIP CVR MD-ASL' , SUBJECT.range_cvr, 'CVR')
% preACZ ATA
SaveDataDICOM(SUBJECT.preACZ.ATA, fullfile(SUBJECT.DICOMdir, SUBJECT.preACZfilenameDCM_ATA), fullfile(SUBJECT.DICOMRESULTSdir, [SUBJECT.preACZfilenameDCM_ATA '.dcm']), 'WIP ATA MD-ASL' , SUBJECT.range_adult_cbf, 'ATA')
% postACZ CBF
SaveDataDICOM(SUBJECT.postACZ.CBF, fullfile(SUBJECT.DICOMdir, SUBJECT.postACZfilenameDCM_CBF), fullfile(SUBJECT.DICOMRESULTSdir, [SUBJECT.postACZfilenameDCM_CBF '.dcm']), 'WIP postACZ CBF MD-ASL', SUBJECT.range_adult_cbf, 'CBF')
% postACZ AAT
SaveDataDICOM(SUBJECT.postACZ.AAT_smth, fullfile(SUBJECT.DICOMdir, SUBJECT.postACZfilenameDCM_AAT), fullfile(SUBJECT.DICOMRESULTSdir, [SUBJECT.postACZfilenameDCM_AAT '.dcm']), 'WIP postACZ AAT MD-ASL', SUBJECT.range_AAT, 'AAT')
% postACZ ATA
SaveDataDICOM(SUBJECT.postACZ.ATA, fullfile(SUBJECT.DICOMdir, SUBJECT.postACZfilenameDCM_ATA), fullfile(SUBJECT.DICOMRESULTSdir, [SUBJECT.postACZfilenameDCM_ATA '.dcm']), 'WIP ATA MD-ASL' , SUBJECT.range_adult_cbf, 'ATA')
disp('CBF, AAT, ATA, CVR Results: DICOMs created')


