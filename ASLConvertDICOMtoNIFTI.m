function ASLConvertDICOMtoNIFTI(DICOMinputdir, NIFTIoutputdir)
% ClinicalASL toolbox 2023, JCWSiero

% % make copy of original DICOM folder (otherwise can't load in VM or scanners)
% system(['mkdir ' DICOMinputdir '/ORIG']);
% 
% % rename DICOM files in DICOM folder, and convert to NIFTI
% system(['dcm2niix -w 0 -r y -f %p_%s ' DICOMinputdir])
% system(['rm -f ' DICOMinputdir '/*_Raw ']);
% system(['rm -f ' DICOMinputdir '/*_PS ']);
% system(['mv -f ' DICOMinputdir '/IM_* ' DICOMinputdir '/ORIG']);
% system(['mv -f ' DICOMinputdir '/XX_* ' DICOMinputdir '/ORIG']);
% system(['mv -f ' DICOMinputdir '/PS_* ' DICOMinputdir '/ORIG']);
% 
% % remove '-' character from filename, otherwise problems with structure fieldnames
% system(['rename -v -f '  '''s/-_//g''' ' ' DICOMinputdir '*']);
% % Convert DICOM files to NIFTI output in NIFTI folder
% system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir]);

% find SOURCE data vTR ASL
files_vTRASL_DCM = dir([DICOMinputdir, '*PRIDE*SOURCE*vTR*']);
if ~isempty(files_vTRASL_DCM)
    system(['rm ' NIFTIoutputdir '/PRIDE*SOURCE*vTR*.*'])
    for i=1:length(files_vTRASL_DCM)
    system(['dcm2niix -w 2 -z y -b y -f %p_%s -s y -o ' NIFTIoutputdir ' ' DICOMinputdir '/' files_vTRASL_DCM(i,1).name]);
    end
end
% remove '-' character from filename, otherwise problems with structure fieldnames
system(['rename -v -f '  '''s/-_//g''' ' ' NIFTIoutputdir '*']);
