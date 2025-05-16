function ASLConvertDICOMtoNIFTI(DICOMinputdir, NIFTIoutputdir, vTR_only)
% ClinicalASL toolbox 2023, JCWSiero

% check operating system: unix vs windows
if isunix  % unix
    % make copy of original DICOM folder (otherwise can't load in VM or scanners)
    system(['mkdir ' fullfile(DICOMinputdir,'ORIG') ' * > /dev/null']);

    % rename DICOM files in DICOM folder, and convert to NIFTI
    system(['dcm2niix -w 0 -r y -f %p_%s ' DICOMinputdir '> dcm2niix_rename.log 2>&1'])
    system(['rm -f ' fullfile(DICOMinputdir,'*._Raw') ' * > /dev/null']);
    system(['rm -f ' fullfile(DICOMinputdir,'*_PS') ' * > /dev/null']);
    system(['mv -f ' fullfile(DICOMinputdir,'IM_*') ' ' fullfile(DICOMinputdir,'ORIG') ' * > /dev/null']);
    system(['mv -f ' fullfile(DICOMinputdir,'XX_*') ' ' fullfile(DICOMinputdir,'ORIG') ' * > /dev/null']);
    system(['mv -f ' fullfile(DICOMinputdir,'PS_*') ' ' fullfile(DICOMinputdir,'ORIG') ' * > /dev/null']);

    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['rename -v -f '  '''s/-_//g''' ' ' fullfile(DICOMinputdir,'*') ' * > /dev/null']);
  
    % Convert DICOM files to NIFTI output in NIFTI folder
    % check for vTR-ASL
    switch nargin
        case 2
            system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir '> dcm2niix.log 2>&1']);
        case 3
            if strcmp(vTR_only,'vTR_only')
                % find SOURCE data vTR ASL
                files_vTRASL_DCM = dir(fullfile(DICOMinputdir, '*PRIDE*SOURCE*vTR*'));
                if ~isempty(files_vTRASL_DCM)
                        system(['rm ' fullfile(NIFTIoutputdir, 'PRIDE*SOURCE*vTR*.*') ' * > /dev/null']);
                    for i=1:length(files_vTRASL_DCM)
                        system(['dcm2niix -w 2 -z y -b y -f %p_%s -s y -o ' NIFTIoutputdir ' ' fullfile(DICOMinputdir, files_vTRASL_DCM(i,1).name) '> dcm2niix.log 2>&1']);
                    end
                end
                 % find SOURCE M0 data vTR ASL
                files_M0_DCM = dir([DICOMinputdir, '*SOURCE*M0*']);
                if ~isempty(files_M0_DCM)
                    system(['del /F /Q ' NIFTIoutputdir '*SOURCE*M0*.*'  ' * > /dev/null']);
                    for i=1:length(files_M0_DCM)
                        system(['dcm2niix -w 2 -z y -b y -f %p_%s -s y -o ' NIFTIoutputdir ' ' DICOMinputdir '\' files_M0_DCM(i,1).name '> dcm2niix.log 2>&1']);
                    end
                end
            end
    end
    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['rename -v -f '  '''s/-_//g''' ' ' fullfile(NIFTIoutputdir,'*') ' * > /dev/null']);
    system(['rm '  fullfile(NIFTIoutputdir,'*Raw*') ' * > /dev/null']);

elseif ispc % windows
    % make copy of original DICOM folder (otherwise can't load in VM or scanners)
    system(['mkdir ' fullfile(DICOMinputdir,'ORIG')]);
 
    % rename DICOM files in DICOM folder, and convert to NIFTI
    system(['dcm2niix -w 0 -r y -f %p_%s ' DICOMinputdir]);
    system(['cmd /c "del /F /Q ' fullfile(DICOMinputdir,'*._Raw') ' "']);
    system(['cmd /c "del /F /Q ' fullfile(DICOMinputdir,'*_PS') ' "']);
    system(['cmd /c "move /Y ' fullfile(DICOMinputdir,'IM_*') ' ' fullfile(DICOMinputdir,'ORIG') ' "']);
    system(['cmd /c "move /Y ' fullfile(DICOMinputdir,'XX_*') ' ' fullfile(DICOMinputdir,'ORIG') ' "']);
    system(['cmd /c "move /Y ' fullfile(DICOMinputdir,'PS_*') ' ' fullfile(DICOMinputdir,'ORIG') ' "']);

    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['cmd /e:on /v:on /c "for %f in ("' DICOMinputdir '\*_-_*") do (set "n=%~nxf" & set "n=!n:_-_=_!" & move /Y "%~ff" ' DICOMinputdir '\"!n!" )"']);
   
    % Convert DICOM files to NIFTI output in NIFTI folder
    % check for vTR-ASL
    switch nargin
        case 2
            system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir '> dcm2niix.log 2>&1']);
        case 3
            if strcmp(vTR_only,'vTR_only')
                % find SOURCE ASL CBF and AAT data vTR ASL
                files_vTRASL_DCM = dir(fullfile(DICOMinputdir, '*PRIDE*SOURCE*vTR*'));
                if ~isempty(files_vTRASL_DCM)
                    system(['del /F /Q ' fullfile(NIFTIoutputdir, 'PRIDE*SOURCE*vTR*.*')]);
                    for i=1:length(files_vTRASL_DCM)
                        system(['dcm2niix -w 2 -z y -b y -f %p_%s -s y -o ' NIFTIoutputdir ' ' fullfile(DICOMinputdir, files_vTRASL_DCM(i,1).name) '> dcm2niix.log 2>&1']);
                    end
                end
                % find SOURCE M0 data vTR ASL
                files_M0_DCM = dir([DICOMinputdir, '*SOURCE*M0*']);
                if ~isempty(files_M0_DCM)
                    system(['del /F /Q ' NIFTIoutputdir '*SOURCE*M0*.*']);
                    for i=1:length(files_M0_DCM)
                        system(['dcm2niix -w 2 -z y -b y -f %p_%s -s y -o ' NIFTIoutputdir ' ' DICOMinputdir '\' files_M0_DCM(i,1).name '> dcm2niix.log 2>&1']);
                    end
                end
            end
    end
    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['cmd /e:on /v:on /c "for %f in ("' NIFTIoutputdir '\*_-_*") do (set "n=%~nxf" & set "n=!n:_-_=_!" & move /Y "%~ff" ' NIFTIoutputdir '\"!n!" )"']);
    system(['del /F /Q ' fullfile(NIFTIoutputdir,'*Raw*')]);

end
