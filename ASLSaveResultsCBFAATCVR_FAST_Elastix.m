function SUBJECT = ASLSaveResultsCBFAATCVR_FAST_Elastix(SUBJECT)
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

elseif strcmp(SUBJECT.RegistrationMethod,'matlab_imreg') % use Matlab imregtform
    optimizer = registration.optimizer.RegularStepGradientDescent;
    metric = registration.metric.MattesMutualInformation;

    disp('Registration T1fromM0 postACZ to preACZ data *********************************************************************')
    fixed_im = SUBJECT.postACZ.M0.*SUBJECT.postACZ.brainmask;
    moving_im = SUBJECT.preACZ.M0.*SUBJECT.preACZ.brainmask;
    tform = imregtform(moving_im, imref3d(size(moving_im)), fixed_im, imref3d(size(fixed_im)), "rigid",optimizer,metric, "PyramidLevels",2);

    moving_reg = imwarp(moving_im, tform, "cubic","OutputView", imref3d(size(moving_im)));
    SaveDataNII(moving_reg, SUBJECT.postACZ_M0_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_M0, scalesclope, [], SUBJECT.TR);
    system(['fslcpgeom ' SUBJECT.preACZ_T1fromM0_path ' ' SUBJECT.postACZ_T1fromM0_2preACZ_path ' -d']);

    disp('Registration CBF postACZ to preACZ *********************************************************************')
    moving_reg = imwarp(SUBJECT.postACZ.CBF, tform, "cubic","OutputView", imref3d(size(SUBJECT.postACZ.CBF)));
    SaveDataNII(moving_reg, SUBJECT.postACZ_CBF_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_CBF, scalesclope, [], SUBJECT.TR);
    system(['fslcpgeom ' SUBJECT.preACZ_CBF_path ' ' postACZ_CBF_2preACZ_path ' -d']);

    disp('Registration AAT postACZ to preACZ *********************************************************************')
    moving_reg = imwarp(SUBJECT.postACZ.AAT, tform, "cubic","OutputView", imref3d(size(SUBJECT.postACZ.AAT)));
    SaveDataNII(moving_reg, postACZ_AAT_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_AAT, scalesclope, [], SUBJECT.TR);
    system(['fslcpgeom ' SUBJECT.preACZ_AAT_path ' ' SUBJECT.postACZ_AAT_2preACZ_path ' -d']);

    disp('Registration ATA postACZ to preACZ *********************************************************************')
    moving_reg = imwarp(SUBJECT.postACZ.ATA, tform, "cubic","OutputView", imref3d(size(SUBJECT.postACZ.ATA)));
    SaveDataNII(moving_reg, postACZ_ATA_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_ATA, scalesclope, [], SUBJECT.TR);
    system(['fslcpgeom ' SUBJECT.preACZ_ATA_path ' ' SUBJECT.postACZ_ATA_2preACZ_path ' -d']);

    disp('Registration mask postACZ to preACZ *********************************************************************')
    moving_reg = imwarp(SUBJECT.postACZ.brainmask, tform, "nearest","OutputView", imref3d(size(SUBJECT.postACZ.CBF)));
    SaveDataNII(moving_reg, SUBJECT.postACZ_mask_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_M0, 1, [], SUBJECT.TR);
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

%% Save to DICOMS: CBF, AAT, ATA, CVR, 
% preACZ CBF
info = dicominfo([SUBJECT.DICOMdir, SUBJECT.preACZfilenameDCM_CBF]); % reference DICOM file
image = SUBJECT.preACZ.CBF;
name = 'WIP preACZ CBF MD-ASL';
dicomname = [SUBJECT.preACZfilenameDCM_CBF '.dcm'];
[a,b,c] = size(image);
scalingfactor = (2^16)/max(image,[],'all'); % for conversion to unsigned int16, divided by 10 otherwize clipping
info.SeriesDescription = name;
info.ProtocolName = name;

for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_adult_cbf);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_adult_cbf);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames = size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.ASLdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.DICOMRESULTSdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% preACZ AAT
info = dicominfo([SUBJECT.DICOMdir, SUBJECT.preACZfilenameDCM_AAT]); % reference DICOM file
image = SUBJECT.preACZ.AAT;
name = 'WIP preACZ AAT(s) MD-ASL';
dicomname = [SUBJECT.preACZfilenameDCM_AAT '.dcm'];
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
info.NumberOfFrames = size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.ASLdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.DICOMRESULTSdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% CVR
info = dicominfo([SUBJECT.DICOMdir, SUBJECT.preACZfilenameDCM_CVR]); % reference DICOM file
image = SUBJECT.CVR_smth;
name = 'WIP CVR MD-ASL';
dicomname = [SUBJECT.preACZfilenameDCM_CVR '.dcm'];
[a,b,c] = size(image);
scalingfactor = (2^15-1)/max(image,[],'all'); % for conversion to signed int16
info.SeriesDescription = name;
info.ProtocolName = name;

for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_cvr);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_cvr);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames = size(image, 3);
dicomwrite(flipud(permute(reshape(int16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.ASLdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);
dicomwrite(flipud(permute(reshape(int16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.DICOMRESULTSdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% preACZ ATA
info = dicominfo([SUBJECT.DICOMdir, SUBJECT.preACZfilenameDCM_ATA]); % reference DICOM file, change when dummy ImageAlgebra Job exists
image = SUBJECT.preACZ.ATA;
name = 'WIP preACZ ATA MD-ASL';
dicomname = [SUBJECT.preACZfilenameDCM_ATA '.dcm'];
[a,b,c] = size(image);
scalingfactor = (2^16)/max(image,[],'all'); % for conversion to unsigned int16, divided by 10 otherwize clipping
info.SeriesDescription = name;
info.ProtocolName = name;
for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_adult_cbf);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_adult_cbf);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames = size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.ASLdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.DICOMRESULTSdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% postACZ CBF
info_ref= dicominfo([SUBJECT.DICOMdir, SUBJECT.postACZfilenameDCM_CBF]); % reference DICOM file
info = info_ref; %copy DICOM info of reference DICOM
image = SUBJECT.postACZ.CBF;
name = 'WIP postACZ CBF MD-ASL';
dicomname = [SUBJECT.postACZfilenameDCM_CBF '.dcm'];
[a,b,c] = size(image);
scalingfactor = (2^16)/max(image,[],'all')/10; % for conversion to unsigned int16, divided by 10 otherwize clipping
info.SeriesDescription = name;
info.ProtocolName = name;
for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_adult_cbf);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_adult_cbf);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames = size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.ASLdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.DICOMRESULTSdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% postACZ AAT
info = dicominfo([SUBJECT.DICOMdir, SUBJECT.postACZfilenameDCM_AAT]); % reference DICOM file
image = SUBJECT.postACZ.AAT;
name = 'WIP postACZ AAT(s) MD-ASL';
dicomname = [SUBJECT.postACZfilenameDCM_AAT '.dcm'];
[a,b,c] = size(image);
scalingfactor = (2^16)/max(image,[],'all')/10; % for conversion to unsigned int16
info.SeriesDescription = name;
info.ProtocolName = name;
for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_AAT);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_AAT);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames = size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.ASLdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.DICOMRESULTSdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% postACZ ATA
info_ref= dicominfo([SUBJECT.DICOMdir, SUBJECT.postACZfilenameDCM_ATA]); % reference DICOM file
info = info_ref; %copy DICOM info of reference DICOM
image = SUBJECT.postACZ.ATA;
name = 'WIP postACZ ATA MD-ASL';
dicomname = [SUBJECT.postACZfilenameDCM_ATA '.dcm'];
[a,b,c] = size(image);
scalingfactor = (2^16)/max(image,[],'all')/10; % for conversion to unsigned int16, divided by 10 otherwize clipping
info.SeriesDescription = name;
info.ProtocolName = name;
for i=1:info.NumberOfFrames
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).PixelValueTransformationSequence.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).Private_2005_140f.Item_1.RescaleSlope = 1/scalingfactor;
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowCenter = mean(SUBJECT.range_adult_cbf);
    info.PerFrameFunctionalGroupsSequence.("Item_"+ num2str(i)).FrameVOILUTSequence.Item_1.WindowWidth = range(SUBJECT.range_adult_cbf);
    if i > size(image, 3) %remove extra frames larger than slice number input image
        info.PerFrameFunctionalGroupsSequence = rmfield(info.PerFrameFunctionalGroupsSequence,("Item_"+ num2str(i)));
    end
end
info.NumberOfFrames = size(image, 3);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.ASLdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.DICOMRESULTSdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

disp('CBF, AAT, ATA, CVR Results: DICOMs created')


