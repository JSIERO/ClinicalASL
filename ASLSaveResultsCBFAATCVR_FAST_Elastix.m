function SUBJECT = ASLSaveResultsCBFAATCVR_FAST_Elastix(SUBJECT)
% ClinicalASL toolbox 2025, JCWSiero
ORprefix='';
%% CBF using 2ndtolast PLD, compute CVR
% register postACZ to preACZ, save BASIl results: NIFTI and PNG - loops over smoothed BASIl data

% Registration postACZ to preACZ space using T1fromM0 AFNI 3dAllineate, 6 dof wsinc interpolation, brain masks as weight
preACZ_T1fromM0_path = [SUBJECT.ASLdir 'preACZ_T1fromM0.nii.gz'];
postACZ_T1fromM0_path = [SUBJECT.ASLdir 'postACZ_T1fromM0.nii.gz'];
postACZ_T1fromM0_2preACZ_path = [SUBJECT.ASLdir 'postACZ_T1fromM0_2preACZ.nii.gz'];
preACZ_mask_path = [SUBJECT.ASLdir 'preACZ_M0_brain_mask.nii.gz'];
postACZ_mask_path = [SUBJECT.ASLdir 'postACZ_M0_brain_mask.nii.gz'];

preACZ_CBF_path = [SUBJECT.ASLdir 'preACZ_BASIL_2tolastPLD_forCBF' ORprefix '/native_space/perfusion_calib.nii.gz'];
postACZ_CBF_path = [SUBJECT.ASLdir 'postACZ_BASIL_2tolastPLD_forCBF' ORprefix '/native_space/perfusion_calib.nii.gz'];
postACZ_CBF_2preACZ_path = [SUBJECT.ASLdir 'postACZ_CBF' ORprefix '_2preACZ.nii.gz'];

preACZ_allPLD_AAT_path = [SUBJECT.ASLdir 'preACZ_BASIL_allPLD_forAAT' ORprefix '/native_space/arrival.nii.gz'];
postACZ_allPLD_AAT_path = [SUBJECT.ASLdir 'postACZ_BASIL_allPLD_forAAT' ORprefix '/native_space/arrival.nii.gz'];
postACZ_allPLD_AAT_2preACZ_path = [SUBJECT.ASLdir 'postACZ_BASIL_allPLD_forAAT' ORprefix '/native_space/arrival_2preACZ.nii.gz'];
postACZ_2preACZ_mask_path = [SUBJECT.ASLdir 'postACZ_M0_brain_mask_2preACZ.nii.gz'];

% Registration postACZ to preACZ
disp('Registration M0 postACZ to preACZ data')
if strcmp(SUBJECT.RegistrationMethod,'elastix') % use elastix
    ElastixParameterFile = fullfile(SUBJECT.GITHUB_ClinicalASLDIR, 'generalFunctions', 'Par0001rigid_6DOF_MI_NIFTIGZ.txt');
    % use Elastix: S. Klein, M. Staring, K. Murphy, M.A. Viergever, J.P.W. Pluim, "elastix: a toolbox for intensity based medical image registration," IEEE Transactions on Medical Imaging, vol. 29, no. 1, pp. 196 - 205, January 2010

    % Registration M0 postACZ to preACZ data
    system(['elastix -f ' preACZ_T1fromM0_path ' -m ' postACZ_T1fromM0_path ' -fMask ' preACZ_mask_path ' -p ' ElastixParameterFile ' -loglevel info -out ' SUBJECT.ASLdir]);
    system(['mv -f ' fullfile(SUBJECT.ASLdir, 'result.0.nii.gz') ' ' postACZ_T1fromM0_2preACZ_path]);
    system(['fslcpgeom ' preACZ_T1fromM0_path ' ' postACZ_T1fromM0_2preACZ_path ' -d']);

    % Registration CBF postACZ to preACZ
    system(['transformix -in ' postACZ_CBF_path ' -out ' SUBJECT.ASLdir ' -tp ' fullfile(SUBJECT.ASLdir, 'TransformParameters.0.txt')]);
    system(['mv -f ' fullfile(SUBJECT.ASLdir, 'result.nii.gz') ' ' postACZ_CBF_2preACZ_path]);
    system(['fslcpgeom ' preACZ_CBF_path ' ' postACZ_CBF_2preACZ_path ' -d']);

    % Registration AAT postACZ to preACZ
    system(['transformix -in ' postACZ_AAT_path ' -out ' SUBJECT.ASLdir ' -tp ' fullfile(SUBJECT.ASLdir, 'TransformParameters.0.txt')]);
    system(['mv -f ' fullfile(SUBJECT.ASLdir, 'result.nii.gz') ' ' postACZ_AAT_2preACZ_path]);
    system(['fslcpgeom ' preACZ_allPLD_AAT_path ' ' postACZ_allPLD_AAT_2preACZ_path ' -d']);

    % Registration mask postACZ to preACZ
    % change ResampleInterpolator to nearestneighbour for binary mask transformation to preACZ
    oldTextLine = '(ResampleInterpolator "FinalBSplineInterpolator")';
    newTextLine = '(ResampleInterpolator "FinalNearestNeighborInterpolator")';
    ChangeElastixParameterFileEntry(fullfile(SUBJECT.ASLdir, 'TransformParameters.0.txt'), oldTextLine, newTextLine, fullfile(SUBJECT.ASLdir, 'TransformParameters.0.NN.txt'));
    system(['transformix -in ' postACZ_mask_path ' -out ' SUBJECT.ASLdir ' -tp ' fullfile(SUBJECT.ASLdir, 'TransformParameters.0.NN.txt')]);
    system(['mv -f ' fullfile(SUBJECT.ASLdir, 'result.nii.gz') ' ' postACZ_mask_2preACZ_path]);
    system(['fslcpgeom ' preACZ_mask_path ' ' postACZ_2preACZ_mask_path ' -d']);
elseif strcmp(SUBJECT.RegistrationMethod,'matlab_imreg') % use Matlab imregtform
    optimizer =  registration.optimizer.RegularStepGradientDescent;
    metric = registration.metric.MattesMutualInformation;

    fixed_im = SUBJECT.postACZ.M0.*SUBJECT.postACZ.brainmask;
    moving_im = SUBJECT.preACZ.M0.*SUBJECT.preACZ.brainmask;
    tform = imregtform(moving_im, imref3d(size(moving_im)), fixed_im, imref3d(size(fixed_im)), "rigid",optimizer,metric, "PyramidLevels",2);

    moving_reg = imwarp(moving_im, tform, "cubic","OutputView", imref3d(size(moving_im)));
    SaveDataNII(moving_reg, postACZ_M0_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_M0, scalesclope, [], SUBJECT.TR);
    system(['fslcpgeom ' preACZ_T1fromM0_path ' ' postACZ_T1fromM0_2preACZ_path ' -d']);

    moving_reg = imwarp(SUBJECT.postACZ.CBF, tform, "cubic","OutputView", imref3d(size(SUBJECT.postACZ.CBF)));
    SaveDataNII(moving_reg, postACZ_CBF_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_CBF, scalesclope, [], SUBJECT.TR);
    system(['fslcpgeom ' preACZ_CBF_path ' ' postACZ_CBF_2preACZ_path ' -d']);

    moving_reg = imwarp(SUBJECT.postACZ.AAT, tform, "cubic","OutputView", imref3d(size(SUBJECT.postACZ.AAT)));
    SaveDataNII(moving_reg, postACZ_AAT_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_AAT, scalesclope, [], SUBJECT.TR);
    system(['fslcpgeom ' preACZ_allPLD_AAT_path ' ' postACZ_allPLD_AAT_2preACZ_path ' -d']);

    moving_reg = imwarp(SUBJECT.postACZ.brainmask, tform, "nearest","OutputView", imref3d(size(SUBJECT.postACZ.CBF)));
    SaveDataNII(moving_reg, postACZ_mask_2preACZ_path(1:end-7), SUBJECT.dummyfilenameSaveNII_M0, 1, [], SUBJECT.TR);
    system(['fslcpgeom ' preACZ_mask_path ' ' postACZ_2preACZ_mask_path ' -d']);
end
disp('Registration finished..')

SlicerPNGs(preACZ_T1fromM0_path, postACZ_T1fromM0_2preACZ_path,  'postACZ_T1fromM0', 'preACZ_T1fromM0', SUBJECT.ASLdir)

% load BASIL CBF results
SUBJECT.preACZ.(['CBF' ORprefix]) = double(niftiread(preACZ_CBF_path));
SUBJECT.postACZ.(['CBF' ORprefix]) = double(niftiread(postACZ_CBF_path));
SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ']) = double(niftiread(postACZ_CBF_2preACZ_path));

% load registered psotaCZ to preACZ brainmask
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
SUBJECT.(['CVR' ORprefix]) = SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ']) - SUBJECT.preACZ.(['CBF' ORprefix]);
SUBJECT.(['CVR' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.(['CVR' ORprefix]).*SUBJECT.nanmask_reg, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE);

%%% save final CBF, CVR NIFTI and .PNGs
smoothloop = {'', '_smth'};
for i=1:length(smoothloop)
    smthprefix = char(smoothloop(i));
    if isempty(smthprefix)
        SaveDataNII(SUBJECT.preACZ.(['CBF' ORprefix smthprefix]), [SUBJECT.ASLdir 'preACZ_CBF' ORprefix smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
        SaveDataNII(SUBJECT.postACZ.(['CBF' ORprefix smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF' ORprefix smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
        SaveDataNII(SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_CBF' ORprefix '_2preACZ' smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);

        SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF' ORprefix smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
        SaveFIGUREtoPNG(SUBJECT.preACZ.(['CBF' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['preACZ_CBF' ORprefix smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
        SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' ORprefix '_2preACZ' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
        SaveFIGUREtoPNG(SUBJECT.postACZ.(['CBF' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['postACZ_CBF' ORprefix '_2preACZ' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
    end
    if strcmp(smthprefix,'_smth')
        SaveDataNII(SUBJECT.(['CVR' ORprefix '_smth']), [SUBJECT.ASLdir 'CVR' ORprefix '_smth'], SUBJECT.dummyfilenameSaveNII,1,[], SUBJECT.TR);
        SaveFIGUREtoPNG(SUBJECT.(['CVR' ORprefix '_smth']), SUBJECT.nanmask_reg, SUBJECT.range_cvr, SUBJECT.RESULTSdir, ['CVR' ORprefix '_smth'],'CVR', 'vik');
    end
end
disp('CBF, CVR Results: NIFTI and .PNGs created')

%% Arrival transit time (AAT), aCBV, arterial transit artefact (ATA) for spatial COV computation, allPLD for another CBF map using all PLD
% register postACZ to preACZ, save BASIl results: NIFTI and PNG
% %%%% % all PLD for AAT (arterial arrival time) map

% load BASIL and registered postACZ to preACZ data
SUBJECT.preACZ.(['AAT' ORprefix]) = double(niftiread(preACZ_allPLD_AAT_path));
SUBJECT.postACZ.(['AAT' ORprefix]) = double(niftiread(postACZ_allPLD_AAT_path));
SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ']) = double(niftiread(postACZ_allPLD_AAT_2preACZ_path));

% Smooth data
SUBJECT.preACZ.(['AAT' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.preACZ.(['AAT' ORprefix]).*SUBJECT.preACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['AAT' ORprefix '_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['AAT' ORprefix]).*SUBJECT.postACZ.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ_smth']) = ASLSmoothImage(SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ']).*SUBJECT.postACZ.nanmask_2preACZ, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

%% save final CBF, CVR NIFTI and .PNGs
smoothloop = {'_smth'};
for i=1:length(smoothloop)
    smthprefix = char(smoothloop(i));
    SaveDataNII(SUBJECT.preACZ.(['AAT' ORprefix smthprefix]), [SUBJECT.ASLdir 'preACZ_AAT' ORprefix smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.postACZ.(['AAT' ORprefix smthprefix]), [SUBJECT.ASLdir 'postACZ_AAT' ORprefix smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ' smthprefix]), [SUBJECT.ASLdir 'postACZ_AAT' ORprefix '_2preACZ' smthprefix], SUBJECT.dummyfilenameSaveNII, 1, [], SUBJECT.TR);

    SaveFIGUREtoPNG(SUBJECT.preACZ.(['AAT' ORprefix smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['preACZ_AAT' ORprefix smthprefix], 'time', 'devon');
    SaveFIGUREtoPNG(SUBJECT.postACZ.(['AAT' ORprefix '_2preACZ' smthprefix]), SUBJECT.nanmask_reg, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['postACZ_AAT' ORprefix '_2preACZ' smthprefix], 'time', 'devon');  
end

%% Save to DICOMS: CBF, AAT, CVR, delta AAT
% preACZ CBF
info = dicominfo([SUBJECT.DICOMdir,SUBJECT.preACZfilenameDCM_CBF]); % reference DICOM file
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

% preACZ AAT
info = dicominfo([SUBJECT.DICOMdir,SUBJECT.preACZfilenameDCM_AAT]); % reference DICOM file
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

% CVR
info = dicominfo([SUBJECT.DICOMdir,SUBJECT.preACZfilenameDCM_CVR]); % reference DICOM file
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

% postACZ CBF
info_ref= dicominfo([SUBJECT.DICOMdir,SUBJECT.postACZfilenameDCM_CBF]); % reference DICOM file
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

% postACZ AAT
info = dicominfo([SUBJECT.DICOMdir,SUBJECT.postACZfilenameDCM_AAT]); % reference DICOM file
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

disp('CBF, AAT, CVR Results: NIFTI, DICOM and .PNGs created')


