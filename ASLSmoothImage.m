function output = ASLSmoothImage(data, spatialdim, FWHM, voxelsize)
% ClinicalASL toolbox 2023, JCWSiero
% can handle NaNs, 2D and 3D gaussian smoothing

sigma = FWHM/2.355;
inplanevoxelsize = voxelsize(1);
data_smooth = zeros(size(data));
if spatialdim == 2
    disp(['Smoothing 2D - can handle NaNs: FWHM(mm) = ' num2str(FWHM)])
    for s=1:size(data, 3)
        if ~isempty(isnan(data))
            filtWidth = 7;
            filtSigma = sigma/inplanevoxelsize;
            imageFilter = fspecial('gaussian',filtWidth,filtSigma);
            data_smooth(:,:,s) = nanconvn(data(:,:,s),imageFilter, 'nanout');
        else
            disp(['Smoothing 2D: FWHM(mm) = ' num2str(FWHM)])
            data_smooth(:,:,s) = imgaussfilt(data(:,:,s),sigma/inplanevoxelsize, 'FilterDomain', 'spatial'); % here also fitler width of 7 is used internally
        end
    end
    
elseif spatialdim == 3
     if ~isempty(isnan(data))
        disp(['Smoothing 3D - can handle NaNs: FWHM(mm) = ' num2str(FWHM)])
        filtWidth = 7;
        filtSigma = sigma/inplanevoxelsize;
        imageFilter = fspecial3('gaussian',filtWidth,filtSigma);
        data_smooth = nanconvn(data,imageFilter, 'nanout');
     else
        disp(['Smoothing 3D: FWHM(mm) = ' num2str(FWHM)])
        data_smooth = imgaussfilt3(data,sigma./voxelsize, 'FilterDomain', 'spatial');        
    end
end

output=double(data_smooth);


