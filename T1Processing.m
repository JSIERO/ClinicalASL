function T1Processing(SUBJECT, filename)
% ClinicalASL toolbox 2023, JCWSiero
if isfile([SUBJECT.ANATOMYdir '/T1_brain_seg_0.nii.gz']) % check if T1 segmentation has been performed
    disp('Segmentation has already been performed');

else % perform T1 anatomy brain extraction and tissue segmentation
    disp('Performing T1 ANATOMY brain extraction and tissue segmentation');
    system(['fslmaths ' SUBJECT.NIFTIdir filename ' ' SUBJECT.ANATOMYdir 'T1.nii.gz']); % move T1 anatomy NIFTI to SUBJECT.ANATOMYdir as T1.nii.gz

    % brain extraction and mask using combination of AFNI and fslmaths
    system(['3dSkullStrip -input ' SUBJECT.ANATOMYdir 'T1.nii.gz' ' -prefix ' SUBJECT.ANATOMYdir 'T1_brain.nii.gz -overwrite']);
    system(['fslmaths ' SUBJECT.ANATOMYdir 'T1_brain.nii.gz -bin -kernel 2D -dilF ' SUBJECT.ANATOMYdir 'T1_brain_mask.nii.gz']);
    system(['fslmaths ' SUBJECT.ANATOMYdir 'T1.nii.gz -mul ' SUBJECT.ANATOMYdir 'T1_brain_mask.nii.gz ' SUBJECT.ANATOMYdir 'T1_brain.nii.gz']);

    % create Slicer PNGs to judge T1 brainmask
    SlicerPNGs([SUBJECT.ANATOMYdir 'T1_brain_mask'], [SUBJECT.ANATOMYdir 'T1'], 'brainmask', 'highres', SUBJECT.ANATOMYdir)

    % FSL FAST to segment tissues and save in previously created locations
    system(['fast -b -g -B -o ' SUBJECT.ANATOMYdir 'T1_brain.nii.gz '  SUBJECT.ANATOMYdir 'T1_brain.nii.gz'])
    system(['fslmaths ' SUBJECT.ANATOMYdir 'T1_brain_restore.nii.gz '  SUBJECT.ANATOMYdir 'T1_brain.nii.gz']); % replace T1_brain with bias field corrected T1_brain_restore

end

