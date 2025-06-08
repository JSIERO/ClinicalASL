function SUBJECT = ASLSaveResultsCBFAAT_vTRASL(SUBJECT)
% ClinicalASL toolbox 2023, JCWSiero

% copy NIFTIs to ASL_vTR folder
% M0
vTR_M0_path = fullfile(SUBJECT.ASLdir, 'vTR_M0.nii.gz');
if isunix
    system(['cp ' SUBJECT.vTRM0filenameNIFTI ' ' vTR_M0_path]);
elseif ispc
    system(['copy /Y ' SUBJECT.vTRM0filenameNIFTI ' ' vTR_M0_path]);
end

info = niftiinfo(vTR_M0_path);
dummy = niftiread(info);
scalesclope = 1;
SUBJECT.vTR.M0 = double(dummy(:,:,:,1))*info.MultiplicativeScaling; % extract first volume
SaveDataNII(SUBJECT.vTR.M0, vTR_M0_path(1:end-7), SUBJECT.dummyfilenameSaveNII_M0, scalesclope, [], SUBJECT.TR);

% make brain masks
vTR_mask_path = fullfile(SUBJECT.ASLdir, 'vTR_M0_brain_mask.nii.gz');

setenv('FSLOUTPUTTYPE','NIFTI_GZ')
system(['bet2 ' vTR_M0_path ' ' vTR_mask_path(1:end-12) ' -m -f 0.5']);

SUBJECT.vTR.brainmask = double(niftiread(vTR_mask_path));

% CBF
vTR_CBF_path = fullfile(SUBJECT.ASLdir, 'vTR_CBF.nii.gz');
if isunix
    system(['cp ' SUBJECT.vTRCBFfilenameNIFTI ' ' vTR_CBF_path]);
elseif ispc
    system(['copy /Y ' SUBJECT.vTRCBFfilenameNIFTI ' ' vTR_CBF_path]);
end

info = niftiinfo(vTR_CBF_path);
dummy = niftiread(info);
scalesclope = 1;
SUBJECT.vTR.CBF = double(dummy(:,:,:,2))*info.MultiplicativeScaling.*SUBJECT.vTR.brainmask; %  extract 2nd volume (corrected CBF, Buxton fit), 1st volume is uncorrected CBF
SaveDataNII(SUBJECT.vTR.CBF, vTR_CBF_path(1:end-7), SUBJECT.dummyfilenameSaveNII_CBF, scalesclope, [], SUBJECT.TR);

% AAT
vTR_AAT_path = fullfile(SUBJECT.ASLdir, 'vTR_AAT.nii.gz');
if isunix
    system(['cp ' SUBJECT.vTRAATfilenameNIFTI ' ' vTR_AAT_path]);
elseif ispc
    system(['copy /Y ' SUBJECT.vTRAATfilenameNIFTI ' ' vTR_AAT_path]);
end

info = niftiinfo(vTR_AAT_path);
dummy = niftiread(info);
scalesclope = 1;
SUBJECT.vTR.AAT = double(dummy(:,:,:,1))./1e3*info.MultiplicativeScaling.*SUBJECT.vTR.brainmask; % convert to s,  extract first volume
SaveDataNII(SUBJECT.vTR.AAT, vTR_AAT_path(1:end-7), SUBJECT.dummyfilenameSaveNII_AAT, scalesclope, [], SUBJECT.TR);

% load CBF results
SUBJECT.vTR.CBF = double(niftiread(vTR_CBF_path));

% load AAT results
SUBJECT.vTR.AAT = double(niftiread(vTR_AAT_path));

% brainmasks

SUBJECT.vTR.nanmask = double(SUBJECT.vTR.brainmask);
SUBJECT.vTR.nanmask(SUBJECT.vTR.nanmask==0) = NaN;

% Smoothing
SUBJECT.vTR.CBF_smth = ASLSmoothImage(SUBJECT.vTR.CBF.*SUBJECT.vTR.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth
SUBJECT.vTR.AAT_smth = ASLSmoothImage(SUBJECT.vTR.AAT.*SUBJECT.vTR.nanmask, 2, SUBJECT.FWHM, SUBJECT.VOXELSIZE); % 2D Smooth

%%% save final CBF, AAT CVR NIFTI and .PNGs
smoothloop = {'', '_smth'};
for i=1:length(smoothloop)
    smthprefix = char(smoothloop(i));
    SaveDataNII(SUBJECT.vTR.(['CBF' smthprefix]), fullfile(SUBJECT.ASLdir, ['vTR_CBF' smthprefix]), SUBJECT.dummyfilenameSaveNII_CBF, 1, [], SUBJECT.TR);
    SaveDataNII(SUBJECT.vTR.(['AAT' smthprefix]), fullfile(SUBJECT.ASLdir, ['vTR_AAT' smthprefix]), SUBJECT.dummyfilenameSaveNII_AAT, 1, [], SUBJECT.TR);

    SaveFIGUREtoPNG(SUBJECT.vTR.(['CBF' smthprefix]), SUBJECT.vTR.nanmask, SUBJECT.range_child_cbf, SUBJECT.RESULTSdir, ['vTR_CBF' smthprefix '_' num2str(SUBJECT.range_child_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.vTR.(['CBF' smthprefix]), SUBJECT.vTR.nanmask, SUBJECT.range_adult_cbf, SUBJECT.RESULTSdir, ['vTR_CBF' smthprefix '_' num2str(SUBJECT.range_adult_cbf(2))], 'CBF', 'viridis');
    SaveFIGUREtoPNG(SUBJECT.vTR.(['AAT' smthprefix]), SUBJECT.vTR.nanmask, SUBJECT.range_AAT, SUBJECT.RESULTSdir, ['vTR_AAT' smthprefix], 'time', 'devon');
end

%% Save to DICOMS: CBF, AAT
% vTR CBF
info = dicominfo(SUBJECT.vTRfilenameDCM_M0); % reference DICOM file
image = SUBJECT.vTR.CBF;
name = 'WIP CBF vTR-ASL';
dicomname = 'vTR_CBF.dcm';
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
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.ASLdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

% vTR AAT
info = dicominfo(SUBJECT.vTRfilenameDCM_M0); % reference DICOM file
image = SUBJECT.vTR.AAT;
name = 'WIP AAT(s) vTR-ASL';
dicomname = 'vTR_AAT.dcm';
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
dicomwrite(flipud(permute(reshape(uint16(image*scalingfactor),[a,b,1,c]),[2,1,3,4])),fullfile(SUBJECT.ASLdir, dicomname),info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

disp('CBF, AAT, Results: NIFTI, DICOM and .PNGs created')
end


