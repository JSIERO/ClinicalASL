function ASLConvertDICOMtoNIFTI(DICOMinputdir, NIFTIoutputdir)
% ClinicalASL toolbox 2023, JCWSiero
% rename DICOM files in DICOM folder, and convert to NIFTI
eval(['!dcm2niix -w 0 -r y -f %p_%s ' DICOMinputdir])
eval(['!rm -f ' DICOMinputdir '/*_Raw' ])
eval(['!rm -f ' DICOMinputdir '/*_PS' ])
eval(['!rm -f ' DICOMinputdir '/IM_*' ])
eval(['!rm -f ' DICOMinputdir '/XX_*' ])
eval(['!rm -f ' DICOMinputdir '/PS_*' ])

% remove '-' character from filename, otherwise problems with structure fieldnames
eval(['!rename -v '  '''s/-_//g''' ' ' DICOMinputdir '*']);
% Convert DICOM files to NIFTI output in NIFTI folder

eval(['!dcm2niix -w 1 -z y -b y -f %p_%s  -o ' NIFTIoutputdir ' ' DICOMinputdir])
% remove '-' character from filename, otherwise problems with structure fieldnames
eval(['!rename -v '  '''s/-_//g''' ' ' NIFTIoutputdir '*']);
