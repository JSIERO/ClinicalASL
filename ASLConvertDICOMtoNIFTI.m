function ASLConvertDICOMtoNIFTI(DICOMinputdir, NIFTIoutputdir)
% ClinicalASL toolbox 2023, JCWSiero
% rename DICOM files in DICOM folder, and convert to NIFTI
eval(['!dcm2niix -w 0 -r y -f %p_%s ' SUBJECT.DICOMdir])
eval(['!rm -f ' SUBJECT.DICOMdir '/*_Raw' ])
eval(['!rm -f ' SUBJECT.DICOMdir '/*_PS' ])
eval(['!rm -f ' SUBJECT.DICOMdir '/IM_*' ])
eval(['!rm -f ' SUBJECT.DICOMdir '/XX_*' ])
eval(['!rm -f ' SUBJECT.DICOMdir '/PS_*' ])

% remove '-' character from filename, otherwise problems with structure fieldnames
eval(['!rename -v '  '''s/-_//g''' ' ' SUBJECT.DICOMdir '*']);
% Convert DICOM files to NIFTI output in NIFTI folder
if numel(dir(SUBJECT.NIFTIdir)) == 2
    eval(['!dcm2niix -w 1 -z y -b y -f %p_%s  -o ' SUBJECT.NIFTIdir ' ' SUBJECT.DICOMdir])
    % remove '-' character from filename, otherwise problems with structure fieldnames
    eval(['!rename -v '  '''s/-_//g''' ' ' SUBJECT.NIFTIdir '*']);
end