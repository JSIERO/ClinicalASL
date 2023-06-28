function SaveDataNII(data, outputfilenm, dummyfilenm, scaleslope, datascale, TR)
% TR is s
[NII, tmpName]=load_untouch_nii(dummyfilenm);
NII.img=data;
if isfloat(data)
    NII.hdr.dime.datatype= 16;
    NII.hdr.dime.bitpix= 32;
end

if strcmp(scaleslope,'samescaling')
    %keep scaling parameters of original(dummy) file
else
    NII.hdr.dime.scl_slope=scaleslope;
    NII.original.hdr.dime.scl_slope=1;
end
if NII.hdr.dime.scl_inter~=0
    disp('WARNING Scaling intercept is not 0, setting to 0 now')
    NII.hdr.dime.scl_inter=0;
end

NII.hdr.dime.dim(5)=size(data,4); % also works for 3D data
NII.original.hdr.dime.dim(5)=size(data,4);
if size(data,4)>1
    NII.hdr.dime.dim(1)=4;
    NII.original.hdr.dime.dim(1)=4;
end

% setting min/max values for viewer (fslview)
if nargin <= 4
    NII.hdr.dime.cal_max=max(data(:));
    NII.hdr.dime.cal_min=min(data(:));
    NII.hdr.dime.glmax=max(data(:));
    NII.hdr.dime.glmin=min(data(:));
else
    if isempty(datascale)
        minscale = min(data(:));
        maxscale = max(data(:));
    else
        minscale = datascale(1);
        maxscale = datascale(2);
    end
    NII.hdr.dime.cal_max=maxscale;
    NII.hdr.dime.cal_min=minscale;
    NII.hdr.dime.glmax=maxscale;
    NII.hdr.dime.glmin=minscale;
    
   % disp(NII.hdr.dime.cal_min)
   % disp(NII.hdr.dime.cal_max)
    
end

if nargin == 6
    NII.hdr.dime.pixdim(5)=TR;% TR is s
    NII.original.hdr.dime.pixdim(5)=TR;
    % also set time_units and space_units to s and mm
    [space_unit, time_unit] = get_units(NII.hdr);
    if space_unit == 1 && time_unit == 1
        % space unit in mm, time_unit in seconds
    elseif space_unit == 1 && time_unit == 1e-3
        NII.hdr.dime.xyzt_units=10; % set xyzt_units to mm and s
    else
        error('XYZT units in NIFTI header are not in mm and s or ms, header check needed!')
    end
end

[p,f,e]=fileparts(outputfilenm);
if strcmp(e, '.nii')
    save_untouch_nii(NII,outputfilenm);
    eval(['!gzip -f ' outputfilenm]); % JCWS gzip .nii
elseif strcmp(e,'.gz')
    outputfilenm=strrep(outputfilenm, '.gz', '');
    save_untouch_nii(NII,outputfilenm);
    eval(['!gzip -f ' outputfilenm]); % JCWS gzip .nii
else
    error('Supply proper extension such as .nii or .nii.gz, ps NIFTI will ALWAYS be gzipped')
end
if isfile([tmpName '.nii']) 
eval(['!rm ' tmpName '.nii'])
end

function [space_unit, time_unit] = get_units(hdr)

switch bitand(hdr.dime.xyzt_units, 7)	% mask with 0x07
    case 1
        space_unit = 1e+3;		% meter, m
    case 3
        space_unit = 1e-3;		% micrometer, um
    otherwise
        space_unit = 1;			% millimeter, mm
end

switch bitand(hdr.dime.xyzt_units, 56)	% mask with 0x38
    case 16
        time_unit = 1e-3;			% millisecond, ms
    case 24
        time_unit = 1e-6;			% microsecond, us
    otherwise
        time_unit = 1;			% second, s
end

return;

