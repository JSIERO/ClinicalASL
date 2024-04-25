function SaveDataNII(data, outputfilenm, dummyfilenm, scaleslope, datarange, TR)
% TR is s

data_info = niftiinfo(dummyfilenm); % read metadata nifti for svaing nifti with same header

if isfloat(data)
    data_info.Datatype = 'double';
    data_info.raw.datatype = 'FLOAT64';    
    data_info.BitsPerPixel = 64;
    data_info.raw.bitpix = 64;
end

if strcmp(scaleslope,'samescaling')
    %keep scaling parameters of original(dummy) file
else
    data_info.MultiplicativeScaling = scaleslope;
    data_info.raw.scl_slope = scaleslope;
end
if data_info.AdditiveOffset~=0
    disp('WARNING Scaling intercept is not 0, setting to 0 now')
    data_info.AdditiveOffset = 0;
    data_info.raw.scl_inter = 0;
end

data_info.ImageSize(4) = size(data,4);
data_info.raw.dim(5) = size(data,4); % also works for 3D data
if size(data,4)>1
    data_info.raw.dim(1) = 4;
end

% setting min/max values for viewer (fslview)
if nargin <= 4
    data_info.DisplayIntensityRange(1) = min(data(:));
    data_info.DisplayIntensityRange(2) = max(data(:));
    data_info.raw.cal_max = max(data(:));
    data_info.raw.cal_min = min(data(:));
else
    if isempty(datarange)
        data_info.DisplayIntensityRange(1) = min(data(:));
        data_info.DisplayIntensityRange(2) = max(data(:));
        data_info.raw.cal_min = min(data(:));
        data_info.raw.cal_max = max(data(:));
    else
        data_info.DisplayIntensityRange(1) = datarange(1);
        data_info.DisplayIntensityRange(2) = datarange(2);
        data_info.raw.cal_min = datarange(1);
        data_info.raw.cal_max = datarange(2);
    end
   
end

if nargin == 6
    data_info.PixelDimensions(4) = TR; %in s
    data_info.raw.pixdim(5) = TR;% in s
    % also set time_units and space_units to s and mm
    if strcmp(data_info.SpaceUnits,'Millimeter') && strcmp(data_info.TimeUnits, 'Second') 
        % space unit in mm, time_unit in seconds
    elseif ~strcmp(data_info.SpaceUnits,'Millimeter') && ~strcmp(data_info.TimeUnits, 'Second') 
        data_info.SpaceUnits = 'Millimeter';
        data_info.TimeUnits =  'Second';
    else
        error('spatial and time units in NIFTI header are not in mm and s or ms, header check needed!')
    end
end

niftiwrite(data,outputfilenm,data_info,'Compressed', true)

end
