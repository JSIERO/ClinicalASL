function SUBJECT = ASLConvertDICOMtoNIFTI(SUBJECT, imager)
% ClinicalASL toolbox 2023, JCWSiero
SUBJECTdir = SUBJECT.SUBJECTdir;
DICOMinputdir = SUBJECT.DICOMdir;
NIFTIoutputdir = SUBJECT.NIFTIdir;

% check operating system: unix vs windows
if isunix  % unix
    if nargin == 2
        if strcmp(imager, 'IMAGER')
            % make copy of original DICOM folder (otherwise can't load in VM or scanners)
            if ~isfolder(fullfile(SUBJECTdir,'DICOMORIG'))
                system(['mkdir ' fullfile(SUBJECTdir,'DICOMORIG')]); % make DICOM folder in outputfolder SUBJECTdir
            end
            system(['cp ' DICOMinputdir '/* ' fullfile(SUBJECTdir,'DICOMORIG')]);  %copy PACS fetched DICOMS to outputfolder SUBJECTdir
            DICOMinputdir = fullfile(SUBJECTdir,'DICOMORIG');% make this the DICOM input folder in output SUBJECTdir for NIFTI conversion
            SUBJECT.DICOMdir = DICOMinputdir; % make fullfile(SUBJECTdir,'DICOM') the new SUBJECT.DICOMdir for extracting scan parameters and such

            if ~isfolder(fullfile(DICOMinputdir,'ORIG'))
                system(['mkdir ' fullfile(DICOMinputdir,'ORIG')]);
            end

            % rename DICOM files in DICOM folder, and convert to NIFTI
            system(['dcm2niix -w 0 -r y -f %p_%s ' DICOMinputdir '> ' NIFTIoutputdir 'dcm2niix_rename.log 2>&1']);
            system(['rm -f ' fullfile(DICOMinputdir,'*._Raw')]);
            system(['rm -f ' fullfile(DICOMinputdir,'*_PS')]);
            system(['mv -f ' fullfile(DICOMinputdir,'IM_*') ' ' fullfile(DICOMinputdir,'ORIG') ' 2>/dev/null']);
            system(['mv -f ' fullfile(DICOMinputdir,'XX_*') ' ' fullfile(DICOMinputdir,'ORIG') ' 2>/dev/null']);
            system(['mv -f ' fullfile(DICOMinputdir,'PS_*') ' ' fullfile(DICOMinputdir,'ORIG') ' 2>/dev/null']);

            % remove '-' character from filename, otherwise problems with structure fieldnames
            system(['rename -v -f '  '''s/-_//g''' ' ' fullfile(DICOMinputdir,'*')]);
            system(['rename -v -f '  '''s/-//g''' ' ' fullfile(DICOMinputdir,'*')]);

            system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir '> ' NIFTIoutputdir 'dcm2niix.log 2>&1']);
        end
    else
        % make copy of original DICOM folder (otherwise can't load in VM or scanners)
        if ~isfolder(fullfile(DICOMinputdir,'ORIG'))
            system(['mkdir ' fullfile(DICOMinputdir,'ORIG')]);
        end

        % rename DICOM files in DICOM folder, and convert to NIFTI
        system(['dcm2niix -w 0 -r y -f %p_%s ' DICOMinputdir '> ' NIFTIoutputdir 'dcm2niix_rename.log 2>&1']);
        system(['rm -f ' fullfile(DICOMinputdir,'*._Raw')]);
        system(['rm -f ' fullfile(DICOMinputdir,'*_PS')]);
        system(['mv -f ' fullfile(DICOMinputdir,'IM_*') ' ' fullfile(DICOMinputdir,'ORIG') ' 2>/dev/null']);
        system(['mv -f ' fullfile(DICOMinputdir,'XX_*') ' ' fullfile(DICOMinputdir,'ORIG') ' 2>/dev/null']);
        system(['mv -f ' fullfile(DICOMinputdir,'PS_*') ' ' fullfile(DICOMinputdir,'ORIG') ' 2>/dev/null']);

        % remove '-' character from filename, otherwise problems with structure fieldnames
        system(['rename -v -f '  '''s/-_//g''' ' ' fullfile(DICOMinputdir,'*')]);
        system(['rename -v -f '  '''s/-//g''' ' ' fullfile(DICOMinputdir,'*')]);

        % Convert DICOM files to NIFTI output in NIFTI folder
        if nargin == 1
            system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir '> ' NIFTIoutputdir 'dcm2niix.log 2>&1']);
        end
    end
    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['rename -v -f '  '''s/-_//g''' ' ' fullfile(NIFTIoutputdir,'*')]);
    system(['rename -v -f '  '''s/-//g''' ' ' fullfile(NIFTIoutputdir,'*')]);    
    system(['rm '  fullfile(NIFTIoutputdir,'*Raw*') ' 2>/dev/null']);

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
    system(['cmd /e:on /v:on /c "for %f in ("' DICOMinputdir '\*-*") do (set "n=%~nxf" & set "n=!n:-=!" & move /Y "%~ff" ' DICOMinputdir '\"!n!" )"']);

    % Convert DICOM files to NIFTI output in NIFTI folder
    if nargin == 2
        system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir '> dcm2niix.log 2>&1']);
    end
    % remove '-' character from filename, otherwise problems with structure fieldnames
    system(['cmd /e:on /v:on /c "for %f in ("' NIFTIoutputdir '\*_-_*") do (set "n=%~nxf" & set "n=!n:_-_=_!" & move /Y "%~ff" ' NIFTIoutputdir '\"!n!" )"']);
    system(['cmd /e:on /v:on /c "for %f in ("' NIFTIoutputdir '\*-*") do (set "n=%~nxf" & set "n=!n:-=!" & move /Y "%~ff" ' NIFTIoutputdir '\"!n!" )"']);
    system(['del /F /Q ' fullfile(NIFTIoutputdir,'*Raw*')]);

end
