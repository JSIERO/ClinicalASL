function ASLConvertDICOMtoNIFTI(DICOMinputdir, NIFTIoutputdir)
% ClinicalASL toolbox 2023, JCWSiero

% make copy of original DICOM folder (otherwise can't load in VM or scanners)
system(['mkdir ' DICOMinputdir '/ORIG']);

% rename DICOM files in DICOM folder, and convert to NIFTI
system(['dcm2niix -w 0 -r y -f %p_%s ' DICOMinputdir])
system(['rm -f ' DICOMinputdir '/*_Raw ']);
system(['rm -f ' DICOMinputdir '/*_PS ']);
system(['mv -f ' DICOMinputdir '/IM_* ' DICOMinputdir '/ORIG']);
system(['mv -f ' DICOMinputdir '/XX_* ' DICOMinputdir '/ORIG']);
system(['mv -f ' DICOMinputdir '/PS_* ' DICOMinputdir '/ORIG']);

% remove '-' character from filename, otherwise problems with structure fieldnames
system(['rename -v -f '  '''s/-_//g''' ' ' DICOMinputdir '*']);
% Convert DICOM files to NIFTI output in NIFTI folder

system(['dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir])
% remove '-' character from filename, otherwise problems with structure fieldnames
system(['rename -v -f '  '''s/-_//g''' ' ' NIFTIoutputdir '*']);
