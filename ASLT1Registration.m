function ASLT1Registration(SUBJECT, prefix)
% ClinicalASL toolbox 2023, JCWSiero

% Register ASL T1fromM0 to T1 anatomy
eval(['! flirt -in ' SUBJECT.ASLdir prefix '_T1fromM0 -out ' SUBJECT.ASLdir prefix '_T1fromM0ASL_2T1 -ref ' SUBJECT.ANATOMYdir 'T1_brain -omat ' SUBJECT.ASLdir prefix '_ASL_2T1.mat  -cost mutualinfo -dof 6 -searchrx -15 15 -searchry -15 15 -searchrz -15 15 -interp trilinear'])

% Inverse transformation matrix -> T1 anatomy to ASL space
eval(['!convert_xfm -omat ' SUBJECT.ASLdir prefix '_ASL_2T1_INV.mat -inverse ' SUBJECT.ASLdir prefix '_ASL_2T1.mat'])

% Transform  T1 segmented tissues to ASL space

% register T1 brain to ASL space
eval(['! flirt -in ' SUBJECT.ANATOMYdir 'T1_brain  -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2T1_INV.mat -out ' SUBJECT.ANATOMYdir 'T1_2ASL_brain_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp sinc'])
% register T1 brain mask to ASL space
eval(['! flirt -in ' SUBJECT.ANATOMYdir 'T1_brain_mask  -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2T1_INV.mat -out ' SUBJECT.ANATOMYdir 'T1_2ASL_brain_mask_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp nearestneighbour']) 
% register T1 CSF mask to ASL space
eval(['! flirt -in ' SUBJECT.ANATOMYdir 'T1_brain_seg_0  -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2T1_INV.mat -out ' SUBJECT.ANATOMYdir 'T1_2ASL_CSF_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp nearestneighbour' ])
% register T1 GM mask to ASL space
eval(['! flirt -in ' SUBJECT.ANATOMYdir 'T1_brain_seg_1  -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2T1_INV.mat -out ' SUBJECT.ANATOMYdir 'T1_2ASL_GM_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp nearestneighbour' ]) 
% register T1 WM mask to ASL space
eval(['! flirt -in ' SUBJECT.ANATOMYdir 'T1_brain_seg_2 -applyxfm -init ' SUBJECT.ASLdir prefix '_ASL_2T1_INV.mat -out ' SUBJECT.ANATOMYdir 'T1_2ASL_WM_' prefix ' -ref ' SUBJECT.ASLdir prefix '_T1fromM0 -paddingsize 0.0 -interp nearestneighbour' ]) 

% Create  Slicer PNGs to judge Registration
SlicerPNGs([SUBJECT.ASLdir prefix '_T1fromM0ASL_2T1'], [SUBJECT.ANATOMYdir 'T1_brain'], 'ASL', 'highres', SUBJECT.ANATOMYdir)

disp(['T1 and tissue segmentation registration to ' prefix ' finished'])
end