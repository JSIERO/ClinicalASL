function ASLMNIRegistration(SUBJECT, prefix)
% ClinicalASL toolbox 2023, JCWSiero
% Register ASL T1fromM0 to MNI 2mm anatomy
system(['flirt -in ' SUBJECT.ASLdir prefix '_T1fromM0 -out ' SUBJECT.ASLdir prefix '_T1fromM0_2MNI -ref ' SUBJECT.MNIdir 'MNI_T1_2mm_brain -omat ' SUBJECT.ASLdir prefix '_ASL_2MNI.mat -cost corratio -dof 12 -interp trilinear']);

% Inverse transformation matrix -> MNI anatomy to ASL space
system(['convert_xfm -omat ' SUBJECT.ASLdir prefix '_ASL_2MNI_INV.mat -inverse ' SUBJECT.ASLdir prefix '_ASL_2MNI.mat']);

% Transform  MNI segmented tissues to ASL space

% register MNI brain to ASL space
system(['flirt -in ' SUBJECT.MNIdir 'MNI_T1_2mm_brain  -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2MNI_INV.mat -out ' SUBJECT.SUBJECTMNIdir 'MNI_2ASL_brain_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp sinc']);
% register MNI brain mask to ASL space
system(['flirt -in ' SUBJECT.MNIdir 'MNI_BRAINMASK_2mm  -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2MNI_INV.mat -out ' SUBJECT.SUBJECTMNIdir 'MNI_2ASL_brain_mask_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp nearestneighbour']); 
% register MNI CSF mask to ASL space
system(['flirt -in ' SUBJECT.MNIdir 'MNI_T1_2mm_brain_seg_0  -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2MNI_INV.mat -out ' SUBJECT.SUBJECTMNIdir 'MNI_2ASL_CSF_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp nearestneighbour' ]);
% register MNI GM mask to ASL space
system(['flirt -in ' SUBJECT.MNIdir 'MNI_T1_2mm_brain_seg_1  -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2MNI_INV.mat -out ' SUBJECT.SUBJECTMNIdir 'MNI_2ASL_GM_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp nearestneighbour' ]); 
% register MNI WM mask to ASL space
system(['flirt -in ' SUBJECT.MNIdir 'MNI_T1_2mm_brain_seg_2 -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2MNI_INV.mat -out ' SUBJECT.SUBJECTMNIdir 'MNI_2ASL_WM_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp nearestneighbour' ]); 

% Create  Slicer PNGs to judge Registration
SlicerPNGs([SUBJECT.ASLdir prefix '_T1fromM0_2MNI'], [SUBJECT.MNIdir 'MNI_T1_2mm_brain'], 'ASL', 'mni', SUBJECT.SUBJECTMNIdir);

disp(['MNI and tissue segmentation registration to ' prefix ' finished'])
end
