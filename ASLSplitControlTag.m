function [CTRL, TAG]=ASLSplitControlTag(CTRLTAG)
% ClinicalASL toolbox 2023, JCWSiero
% interleave CTRL and TAG data; CTRL1..TAG1..CTRL2..TAG2..
CTRL=[];
TAG=[];
dims = size(CTRLTAG);% xdim y xim zdim tdim asl_label

odd_index=2*(1:dims(4)/2) - 1;
even_index=2*(1:dims(4)/2);

CTRL = CTRLTAG(:,:,:,odd_index);
TAG = CTRLTAG(:,:,:,even_index);

    

