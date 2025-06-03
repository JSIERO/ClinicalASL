function SaveDataDICOM(image, template_dicom_path, output_dicom_path, name, value_range, content_label)
% Export ASL-derived image as multiframe DICOM using Philips reference
%
% Parameters:
%   image                - 3D matrix of ASL data (e.g. CBF, CVR)
%   template_dicom_path  - path to a template Philips multiframe DICOM
%   output_dicom_path    - path to output DICOM
%   name                 - series/protocol name
%   value_range          - [min max] for VOI LUT (windowing)
%   content_label        - e.g. 'CBF', 'CVR', 'AAT' (used to tag DICOM as QUANTITATIVE)

if nargin < 6 || isempty(content_label)
    warning('No content_label specified. Please supply one of: ''CBF'', ''CVR'', ''AAT'', or ''ATA'' to ensure correct DICOM labeling and scaling behavior.');
    return;
end

% Determine units and int type
use_signed = strcmpi(content_label, 'CVR');
switch upper(content_label)
    case {'CBF', 'CVR', 'ATA'}
        unit_str = 'ml/100g/min';
    case 'AAT'
        unit_str = 's';
    otherwise
        unit_str = '';
end

info = dicominfo(template_dicom_path);
[a, b, c] = size(image);

% Scaling and cast
if use_signed
    abs_max = max(abs(image(:)));
    scalingfactor = (2^15 - 1) / abs_max;
    image_scaled = int16(image * scalingfactor);
    pixel_representation = 1;
else
    scalingfactor = (2^16 - 1) / max(image(:));
    image_scaled = uint16(image * scalingfactor);
    pixel_representation = 0;
end

% Update DICOM metadata
info.SeriesInstanceUID     = dicomuid;
info.SOPClassUID           = '1.2.840.10008.5.1.4.1.1.4.1';
info.SeriesDescription     = [name ' [' unit_str ']'];
info.ProtocolName          = name;
info.ContentDescription    = ['ASL derived map: ' upper(content_label) ' (' unit_str ')'];
info.ContentLabel          = upper(content_label);
info.ImageType             = {'DERIVED', 'SECONDARY', 'QUANTITATIVE'};
info.NumberOfFrames        = c;
info.RescaleSlope          = 1 / scalingfactor;
info.RescaleIntercept      = 0;
info.PixelRepresentation   = pixel_representation;
info.BitsAllocated         = 16;
info.BitsStored            = 16;
info.HighBit               = 15;

% Automatically create ReferencedSeriesSequence from reference
info.ReferencedSeriesSequence = struct('ReferencedSeriesInstanceUID', info.SeriesInstanceUID, 'ReferencedInstanceSequence', struct([]));

% Per-frame updates
for i = 1:c
    frame = info.PerFrameFunctionalGroupsSequence.("Item_" + i);
    frame.PixelValueTransformationSequence.Item_1.RescaleSlope = 1 / scalingfactor;
    frame.Private_2005_140f.Item_1.RescaleSlope                = 1 / scalingfactor;
    frame.FrameVOILUTSequence.Item_1.WindowCenter              = mean(value_range);
    frame.FrameVOILUTSequence.Item_1.WindowWidth               = diff(value_range);
    info.PerFrameFunctionalGroupsSequence.("Item_" + i) = frame;
end

% Orient and write
pixeldata = flipud(permute(reshape(image_scaled, [a, b, 1, c]), [2, 1, 3, 4])); % SANINTY CHECK THIS!
dicomwrite(pixeldata, output_dicom_path, info, 'CreateMode', 'Copy', 'MultiframeSingleFile', true);

fprintf('✅ DICOM written to: %s\n  Content: %s (%s), Slope: %.6f\n', output_dicom_path, upper(content_label), unit_str, 1 / scalingfactor);
end
