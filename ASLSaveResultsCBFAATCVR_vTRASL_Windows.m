function SUBJECT = ASLSaveResultsCBFAATCVR_vTRASL_Windows(SUBJECT)
% ClinicalASL toolbox 2023, JCWSiero

% copy NIFTIs to ASL_vTR folder
% M0
preACZ_M0_path = [SUBJECT.ASLdir 'preACZ_M0.nii.gz'];
postACZ_M0_path = [SUBJECT.ASLdir 'postACZ_M0.nii.gz'];
if isunix
    system(['cp ' SUBJECT.preACZM0filenameNIFTI ' ' preACZ_M0_path]);
    system(['cp ' SUBJECT.postACZM0filenameNIFTI ' ' postACZ_M0_path]);
elseif ispc
    system(['copy /Y ' SUBJECT.preACZM0filenameNIFTI ' ' preACZ_M0_path]);
    system(['copy /Y ' SUBJECT.postACZM0filenameNIFTI ' ' postACZ_M0_path]);
end

info = niftiinfo(preACZ_M0_path);
dummy = niftiread(info);
scalesclope = 1;
SUBJECT.preACZ.M0 = double(dummy(:,:,:,1))*info.MultiplicativeScaling; % extract first volume
SaveDataNII(SUBJECT.preACZ.M0, preACZ_M0_path(1:end-7), SUBJECT.dummyfilenameSaveNII_M0, scalesclope, [], SUBJECT.TR);

info = niftiinfo(postACZ_M0_path);
dummy = niftiread(info);
scalesclope = 1;
SUBJECT.postACZ.M0 = double(dummy(:,:,:,1))*info.MultiplicativeScaling;
SaveDataNII(SUBJECT.postACZ.M0, postACZ_M0_path(1:end-7), SUBJECT.dummyfilenameSaveNII_M0, scalesclope, [], SUBJECT.TR);

% make brain masks
preACZ_mask_path = [SUBJECT.ASLdir 'preACZ_M0_brain_mask.nii.gz'];
postACZ_mask_path = [SUBJECT.ASLdir 'postACZ_M0_brain_mask.nii.gz'];

if isunix
    disp('maken!')
elseif ispc
    setenv('FSLOUTPUTTYPE','NIFTI_GZ')
    system(['bet2 ' preACZ_M0_path ' ' preACZ_mask_path(1:end-12) ' -m -f 0.5']);
    system(['bet2 ' postACZ_M0_path ' ' postACZ_mask_path(1:end-12) ' -m -f 0.5']);
end

SUBJECT.preACZ.brainmask = double(niftiread(preACZ_mask_path));
SUBJECT.postACZ.brainmask = double(niftiread(postACZ_mask_path));

% CBF
preACZ_CBF_path = [SUBJECT.ASLdir 'preACZ_CBF.nii.gz'];
postACZ_CBF_path = [SUBJECT.ASLdir 'postACZ_CBF.nii.gz'];
if isunix
    system(['cp ' SUBJECT.preACZCBFfilenameNIFTI ' ' preACZ_CBF_path]);
    system(['cp ' SUBJECT.postACZCBFfilenameNIFTI ' ' postACZ_CBF_path]);
elseif ispc
    system(['copy /Y ' SUBJECT.preACZCBFfilenameNIFTI ' ' preACZ_CBF_path]);
    system(['copy /Y ' SUBJECT.postACZCBFfilenameNIFTI ' ' postACZ_CBF_path]);
end

info = niftiinfo(preACZ_CBF_path);
dummy = niftiread(info);
scalesclope = 1;
SUBJECT.preACZ.CBF = double(dummy(:,:,:,2))*info.MultiplicativeScaling.*SUBJECT.preACZ.brainmask; %  extract 2nd volume (corrected CBF, Buxton fit), 1st volume is uncorrected CBF
SaveDataNII(SUBJECT.preACZ.CBF, preACZ_CBF_path(1:end-7), SUBJECT.dummyfilenameSaveNII_CBF, scalesclope, [], SUBJECT.TR);

info = niftiinfo(postACZ_CBF_path);
dummy = niftiread(info);
scalesclope = 1;
SUBJECT.postACZ.CBF= double(dummy(:,:,:,2))*info.MultiplicativeScaling.*SUBJECT.postACZ.brainmask;
SaveDataNII(SUBJECT.postACZ.CBF, postACZ_CBF_path(1:end-7), SUBJECT.dummyfilenameSaveNII_CBF, scalesclope, [], SUBJECT.TR);

% AAT
preACZ_AAT_path = [SUBJECT.ASLdir 'preACZ_AAT.nii.gz'];
postACZ_AAT_path = [SUBJECT.ASLdir 'postACZ_AAT.nii.gz'];
if isunix
    system(['cp ' SUBJECT.preACZAATfilenameNIFTI ' ' preACZ_AAT_path]);
    system(['cp ' SUBJECT.postACZAATfilenameNIFTI ' ' postACZ_AAT_path]);
elseif ispc
    system(['copy /Y ' SUBJECT.preACZAATfilenameNIFTI ' ' preACZ_AAT_path]);
    system(['copy /Y ' SUBJECT.postACZAATfilenameNIFTI ' ' postACZ_AAT_path]);
end

info = niftiinfo(preACZ_AAT_path);
dummy = niftiread(info);
scalesclope = 1;
SUBJECT.preACZ.AAT = double(dummy(:,:,:,1))./1e3*info.MultiplicativeScaling.*SUBJECT.preACZ.brainmask; % convert to s,  extract first volume
SaveDataNII(SUBJECT.preACZ.AAT, preACZ_AAT_path(1:end-7), SUBJECT.dummyfilenameSaveNII_AAT, scalesclope, [], SUBJECT.TR);

info = niftiinfo(postACZ_AAT_path);
dummy = niftiread(info);
scalesclope = 1;
SUBJECT.postACZ.AAT = double(dummy(:,:,:,1))./1e3*info.MultiplicativeScaling.*SUBJECT.postACZ.brainmask; % convert to s,  extract first volume
SaveDataNII(SUBJECT.postACZ.AAT, postACZ_AAT_path(1:end-7), SUBJECT.dummyfilenameSaveNII_AAT, scalesclope, [], SUBJECT.TR);


%% compute CVR
% Registration postACZ to preACZ space using M0 

postACZ_M0_2preACZ_path = [SUBJECT.ASLdir 'postACZ_M0_2preACZ.nii.gz'];
postACZ_CBF_2preACZ_path = [SUBJECT.ASLdir 'postACZ_CBF_2preACZ.nii.gz'];
postACZ_AAT_2preACZ_path = [SUBJECT.ASLdir 'postACZ_AAT_2preACZ.nii.gz'];
postACZ_mask_2preACZ_path = [SUBJECT.ASLdir 'postACZ_M0_brain_mask_2preACZ.nii.gz'];

disp('Registration M0 postACZ to preACZ data')
if strcmp(SUBJECT.RegistrationMethod,'elastix') % use elastix
    ElastixParameterFile = [SUBJECT.GITHUB_ClinicalASLDIR '\generalFunctions\' 'Par0001rigid_6DOF_MI_NIFTIGZ.txt'];
    % use Elastix: S. Klein, M. Staring, K. Murphy, M.A. Viergever, J.P.W. Pluim, "elastix: a toolbox for intensity based medical image registration," IEEE Transactions on Medical Imaging, vol. 29, no. 1, pp. 196 - 205, January 2010

    % Registration M0 postACZ to preACZ data
    system(['elastix -f ' preACZ_M0_path ' -m ' postACZ_M0_path ' -fMask ' preACZ_mask_path ' -p ' ElastixParameterFile ' -loglevel info -out ' SUBJECT.ASLdir]);
    system(['move /Y ' SUBJECT.ASLdir 'result.0.nii.gz ' postACZ_M0_2preACZ_path]);

    % Registration CBF postACZ to preACZ
    system(['transformix -in ' postACZ_CBF_path ' -out ' SUBJECT.ASLdir ' -tp ' SUBJECT.ASLdir 'TransformParameters.0.txt']);
    system(['move /Y ' SUBJECT.ASLdir 'result.nii.gz ' postACZ_CBF_2preACZ_path]);

    % Registration AAT postACZ to preACZ
    system(['transformix -in ' postACZ_AAT_path ' -out ' SUBJECT.ASLdir ' -tp ' SUBJECT.ASLdir 'TransformParameters.0.txt']);
    system(['move /Y ' SUBJECT.ASLdir 'result.nii.gz ' postACZ_AAT_2preACZ_path]);

    % Registration mask postACZ to preACZ
    % change ResampleInterpolator to nearestneighbour for binary mask transformation to preACZ
    oldTextLine = '(ResampleInterpolator "FinalBSplineInterpolator")';
    newTextLine = '(ResampleInterpolator "FinalNearestNeighborInterpolator")';
    ChangeElastixParameterFileEntry([SUBJECT.ASLdir 'TransformParameters.0.txt'], oldTextLine, newTextLine, [SUBJECT.ASLdir 'TransformParameters.0.NN.txt']);

    system(['transformix -in ' postACZ_mask_path ' -out ' SUBJECT.ASLdir ' -tp ' SUBJECT.ASLdir 'TransformParameters.0.NN.txt']);
    system(['move /Y ' SUBJECT.ASLdir 'result.nii.gz ' postACZ_mask_2preACZ_path]);

elseif strcmp(SUBJECT.RegistrationMethod,'matlab_imreg') % use Matlab imregtform
    optimizer =  registration.optimizer.RegularStepGradientDescent;
    metric = registration.metric.MattesMutualInformation;

    fixed_im = SUBJECT.postACZ.M0.*SUBJECT.postACZ.brainmask;
    moving_im = SUBJECT.preACZ.M0.*SUBJECT.preACZ.brainmask;
    tform = imregtform(moving_im, imref3d(size(moving_im)), fixed_im, imref3d(size(fixed_im)), "rigid",optimizer,metric, "PyramidLevels",2);

    moving_reg = imwarp(moving_im, tform, "cubic","OutputView", imref3d(size(moving_im)));
    SaveDataNII(moving_reg, postACZ_M0_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_M0, scalesclope, [], SUBJECT.TR);

    moving_reg = imwarp(SUBJECT.postACZ.CBF, tform, "cubic","OutputView", imref3d(size(SUBJECT.postACZ.CBF)));
    SaveDataNII(moving_reg, postACZ_CBF_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_CBF, scalesclope, [], SUBJECT.TR);

    moving_reg = imwarp(SUBJECT.postACZ.AAT, tform, "cubic","OutputView", imref3d(size(SUBJECT.postACZ.AAT)));
    SaveDataNII(moving_reg, postACZ_AAT_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_AAT, scalesclope, [], SUBJECT.TR);

    moving_reg = imwarp(SUBJECT.postACZ.brainmask, tform, "nearest","OutputView", imref3d(size(SUBJECT.postACZ.CBF)));
    SaveDataNII(moving_reg, postACZ_mask_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_M0, 1, [], SUBJECT.TR);
end
disp('Registration finished..')

% load CBF results
SUBJECT.preACZ.CBF = double(niftiread(preACZ_CBF_path));
SUBJECT.postACZ.CBF = double(niftiread(postACZ_CBF_path));
SUBJECT.postACZ.CBF_2preACZ = double(niftiread(postACZ_CBF_2preACZ_path));

% load AAT results
SUBJECT.preACZ.AAT = double(niftiread(preACZ_AAT_path));
SUBJECT.postACZ.AAT = double(niftiread(postACZ_AAT_path));
SUBJECT.postACZ.AAT_2preACZ = double(niftiread(postACZ_AAT_2preACZ_path));

% brainmasks
SUBJECT.postACZ.brainmask_2preACZ = double(niftiread(postACZ_mask_2preACZ_path));

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

% Smoothing
SUBJECT.preACZ.CBF_smth = ASLSmoothImage(SUBJECT.preACZ.CBF.*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.CBF_smth = ASLSmoothImage(SUBJECT.postACZ.CBF.*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.CBF_2preACZ_smth = ASLSmoothImage(SUBJECT.postACZ.CBF_2preACZ.*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.preACZ.AAT_smth = ASLSmoothImage(SUBJECT.preACZ.AAT.*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.AAT_smth = ASLSmoothImage(SUBJECT.postACZ.AAT.*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.AAT_2preACZ_smth = ASLSmoothImage(SUBJECT.postACZ.AAT_2preACZ.*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

SUBJECT.CVR_smth = ASLSmoothImage(SUBJECT.CVR.*SUBJECT.nanmask_reg, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE);
% Compute CVR in percentage
SUBJECT.CVR_percentage_smth = SUBJECT.CVR_smth ./SUBJECT.preACZ.CBF_smth .*SUBJECT.nanmask_reg * 100;

% compute deltaAAT map to reflect changes in AAT
SUBJECT.AATdelta_smth = ASLSmoothImage((SUBJECT.postACZ.AAT_2preACZ_smth - SUBJECT.preACZ.AAT_smth).*SUBJECT.nanmask_reg, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE);

%%% save final CBF, AAT CVR NIFTI and .PNGs
smoothloop = {'', '_smth'};
for i=1:length(smoothloop)
    smthprefix = char(smoothloop(i));
    SaveDataNII(SUBJECT.preACZ.(['CBF' smthprefix]), [SUBJECT.ASLdir 'preACZ_CBF' smthprefix], SUBJECT.dummyfilenameSaveNII_CBF, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.postACZ.(['CBF' smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF' smthprefix], SUBJECT.dummyfilenameSaveNII_CBF, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.postACZ.(['CBF_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF_2preACZ' smthprefix], SUBJECT.dummyfilenameSaveNII_CBF, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.preACZ.(['AAT' smthprefix]), [SUBJECT.ASLdir 'preACZ_AAT' smthprefix], SUBJECT.dummyfilenameSaveNII_AAT, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.postACZ.(['AAT' smthprefix]), [SUBJECT.ASLdir 'postACZ_AAT' smthprefix], SUBJECT.dummyfilenameSaveNII_AAT, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.postACZ.(['AAT_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_AAT_2preACZ' smthprefix], SUBJECT.dummyfilenameSaveNII_AAT, 1, [], SUBJECT.TR);

    SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF_2preACZ' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF_2preACZ' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.preACZ.(['AAT' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['preACZ_AAT' smthprefix], 'time', 'devon');
    SaveFIGUREtoPNG(SUBJECT.postACZ.(['AAT' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['postACZ_AAT' smthprefix], 'time', 'devon');
    SaveFIGUREtoPNG(SUBJECT.postACZ.(['AAT_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['postACZ_AAT_2preACZ' smthprefix], 'time', 'devon');
    if strcmp(smthprefix,'_smth')
        SaveDataNII(SUBJECT.CVR_smth, [SUBJECT.ASLdir 'CVR_smth'], SUBJECT.dummyfilenameSaveNII_CBF,1,[], SUBJECT.TR);
        SaveDataNII(SUBJECT.CVR_percentage_smth,[SUBJECT.ASLdir 'CVR_percentage_smth'], SUBJECT.dummyfilenameSaveNII_CBF,1,[], SUBJECT.TR);
        SaveDataNII(SUBJECT.AATdelta_smth, [SUBJECT.ASLdir 'AATdelta_smth'], SUBJECT.dummyfilenameSaveNII_AAT,1,[], SUBJECT.TR);

        SaveFIGUREtoPNG(SUBJECT.CVR_smth, SUBJECT.nanmask_reg, SUBJECT.range_cvr, SUBJECT.RESULTSdir, 'CVR_smth','CVR', 'vik');
        SaveFIGUREtoPNG(SUBJECT.CVR_percentage_smth, SUBJECT.nanmask_reg, [-100 100], SUBJECT.RESULTSdir, 'CVR_PERCENTAGE_smth','%','vik');
        SaveFIGUREtoPNG(SUBJECT.AATdelta_smth, SUBJECT.nanmask_reg, SUBJECT.range_AATdelta, SUBJECT.RESULTSdir, 'AATdelta_smth','time_delta', 'vik');
    end
end

%% Save to DICOMS: CBF, AAT, CVR, delta AAT
% CVR
info = dicominfo(SUBJECT.preACZfilenameDCM_M0); % reference DICOM file
image = SUBJECT.CVR_smth;
name = 'WIP CVR vTR-ASL';
dicomname = 'CVR_smth.dcm';
[a,b,c] = size(image);
scalingfactor = (2^15-1)/max(image,[],'all'); % for conversion to signed int16
info.SeriesDescription = name;
info.ProtocolName = name;
info.SeriesNumber = info.SeriesNumber + 1;

for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_cvr);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_cvr);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames= size(image, 3);
dicomwrite(flipud(permute(reshape(int16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),[SUBJECT.ASLdir dicomname],info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% preACZ CBF
info = dicominfo(SUBJECT.preACZfilenameDCM_M0); % reference DICOM file
image = SUBJECT.preACZ.CBF;
name = 'WIP preACZ CBF vTR-ASL';
dicomname = 'preACZ_CBF.dcm';
[a,b,c] = size(image);
scalingfactor = (2^16)/max(image,[],'all'); % for conversion to unsigned int16, divided by 10 otherwize clipping
info.SeriesDescription = name;
info.ProtocolName = name;
info.SeriesNumber = info.SeriesNumber + 2;

for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_adult_cbf);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_adult_cbf);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames= size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),[SUBJECT.ASLdir dicomname],info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% postACZ CBF
info = dicominfo(SUBJECT.preACZfilenameDCM_M0); % reference DICOM file
image = SUBJECT.postACZ.CBF;
name = 'WIP postACZ CBF vTR-ASL';
dicomname = 'postACZ_CBF.dcm';
[a,b,c] = size(image);
scalingfactor = (2^16)/max(image,[],'all')/10; % for conversion to unsigned int16, divided by 10 otherwize clipping
info.SeriesDescription = name;
info.ProtocolName = name;
info.SeriesNumber = info.SeriesNumber + 3;
for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_adult_cbf);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_adult_cbf);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames= size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),[SUBJECT.ASLdir dicomname],info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% preACZ AAT
info = dicominfo(SUBJECT.preACZfilenameDCM_M0); % reference DICOM file
image = SUBJECT.preACZ.AAT;
name = 'WIP preACZ AAT(s) vTR-ASL';
dicomname = 'preACZ_AAT.dcm';
info.SeriesNumber = info.SeriesNumber + 4;
[a,b,c] = size(image);
scalingfactor = (2^16)/max(image,[],'all')/10; % for conversion to unsigned int16, divided by 10 otherwize clipping
info.SeriesDescription = name;
info.ProtocolName = name;
for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_AAT);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_AAT);
    if i > c %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames= size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),[SUBJECT.ASLdir dicomname],info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% postACZ AAT
info = dicominfo(SUBJECT.preACZfilenameDCM_M0); % reference DICOM file
image = SUBJECT.postACZ.AAT;
name = 'WIP postACZ AAT(s) vTR-ASL';
dicomname = 'postACZ_AAT.dcm';
[a,b,c] = size(image);
scalingfactor = (2^16)/max(image,[],'all')/10; % for conversion to unsigned int16
info.SeriesDescription = name;
info.ProtocolName = name;
info.SeriesNumber = info.SeriesNumber + 5;
for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_AAT);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_AAT);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames= size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),[SUBJECT.ASLdir dicomname],info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% deltaAAT
info = dicominfo(SUBJECT.preACZfilenameDCM_M0); % reference DICOM file
image = SUBJECT.AATdelta_smth;
name = 'WIP deltaAAT(s) vTR-ASL';
dicomname = 'AATdelta_smth.dcm';
[a,b,c] = size(image);
scalingfactor = (2^15-1)/max(image,[],'all'); % for conversion to signed int16
info.SeriesDescription = name;
info.ProtocolName = name;
info.SeriesNumber = info.SeriesNumber + 6;
for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_AATdelta);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_AATdelta);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames= size(image, 3);
dicomwrite(flipud(permute(reshape(int16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),[SUBJECT.ASLdir dicomname],info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

disp('CBF, AAT, CVR Results: NIFTI, DICOM and .PNGs created')
end


