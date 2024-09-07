function ASLConvertDICOMtoNIFTI(DICOMinputdir, NIFTIoutputdir, vTR_only)
% ClinicalASL toolbox 2023, JCWSiero

% check operating system: unix vs windows
if isunix  % unix
    % make copy of original DICOM folder (otherwise can't load in VM or scanners)
    system(['mkdir ' DICOMinputdir 'ORIG']);

    % rename DICOM files in DICOM folder, and convert to NIFTI
    system(['dcm2niix -w 0 -r y -f %p_%s ' DICOMinputdir])
    system(['rm -f ' DICOMinputdir '*_Raw ']);
    system(['rm -f ' DICOMinputdir '*_PS ']);
    system(['mv -f ' DICOMinputdir 'IM_* ' DICOMinputdir 'ORIG']);
    system(['mv -f ' DICOMinputdir 'XX_* ' DICOMinputdir 'ORIG']);
    system(['mv -f ' DICOMinputdir 'PS_* ' DICOMinputdir 'ORIG']);

    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['rename -v -f '  '''s/-_//g''' ' ' DICOMinputdir '*']);
    % Convert DICOM files to NIFTI output in NIFTI folder
    system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir]);
      
    % Convert DICOM files to NIFTI output in NIFTI folder
    % check for vTR-ASL
    switch nargin
        case 2
            system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir]);
        case 3
            if strcmp(vTR_only,'vTR_only')
                % find SOURCE data vTR ASL
                files_vTRASL_DCM = dir([DICOMinputdir, '*PRIDE*SOURCE*vTR*']);
                if ~isempty(files_vTRASL_DCM)
                        system(['rm ' NIFTIoutputdir 'PRIDE*SOURCE*vTR*.*']);
                    for i=1:length(files_vTRASL_DCM)
                        system(['dcm2niix -w 2 -z y -b y -f %p_%s -s y -o ' NIFTIoutputdir ' ' DICOMinputdir '\' files_vTRASL_DCM(i,1).name]);
                    end
                end
            end
    end
    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['rename -v -f '  '''s/-_//g''' ' ' NIFTIoutputdir '*']);

elseif ispc % windows
    % make copy of original DICOM folder (otherwise can't load in VM or scanners)
    system(['mkdir ' DICOMinputdir 'ORIG']);
 
    % rename DICOM files in DICOM folder, and convert to NIFTI
    system(['dcm2niix -w 0 -r y -f %p_%s ' DICOMinputdir]);
    system(['cmd /c "del /F /Q ' DICOMinputdir '*_Raw "']);
    system(['cmd /c "del /F /Q ' DICOMinputdir '*_PS "']);
    system(['cmd /c "move /Y ' DICOMinputdir 'IM_* ' DICOMinputdir 'ORIG"']);
    system(['cmd /c "move /Y ' DICOMinputdir 'XX_* ' DICOMinputdir 'ORIG"']);
    system(['cmd /c "move /Y ' DICOMinputdir 'PS_* ' DICOMinputdir 'ORIG"']);

    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['cmd /e:on /v:on /c "for %f in ("' DICOMinputdir '*_-_*") do (set "n=%~nxf" & set "n=!n:_-_=_!" & move /Y "%~ff" ' DICOMinputdir '"!n!" )"']);
   
    % Convert DICOM files to NIFTI output in NIFTI folder
    % check for vTR-ASL
    switch nargin
        case 2
            system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir]);
        case 3
            if strcmp(vTR_only,'vTR_only')
                % find SOURCE data vTR ASL
                files_vTRASL_DCM = dir([DICOMinputdir, '*PRIDE*SOURCE*vTR*']);
                if ~isempty(files_vTRASL_DCM)
                    system(['del /F /Q ' NIFTIoutputdir 'PRIDE*SOURCE*vTR*.*']);
                    for i=1:length(files_vTRASL_DCM)
                        system(['dcm2niix -w 2 -z y -b y -f %p_%s -s y -o ' NIFTIoutputdir ' ' DICOMinputdir '\' files_vTRASL_DCM(i,1).name]);
                    end
                end
            end
    end
    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['cmd /e:on /v:on /c "for %f in ("' NIFTIoutputdir '*_-_*") do (set "n=%~nxf" & set "n=!n:_-_=_!" & move /Y "%~ff" ' NIFTIoutputdir '"!n!" )"']);
end
