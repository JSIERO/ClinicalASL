 function ASLBASILanalysis(SUBJECT, locationASLlabelcontrolPLDNIFTI, locationM0, locationMask, outputmap, PLDlist, locationBASILinfo, artoff, spatialoff)
% ClinicalASL toolbox 2023, JCWSiero
PLDstring = sprintf('%.05g,' , SUBJECT.PLDS(PLDlist));
PLDstring = PLDstring(1:end-1);% strip final comma

if nargin >= 7
    if isempty(numel(locationBASILinfo))
        locationBASILinfostring=[];
    else
        locationBASILinfostring=[' --model-options=' locationBASILinfo];
    end
else
    locationBASILinfostring=[];
end

if nargin >= 8
    if strcmp (artoff, 'artoff')
        artoffstring= ' --artoff'; % do not inver arterial aCBV component, needed for aterial transit time aterfact (ATA) mapping - also much faster, use for per dynamic CBF analysis (outlierremoval)
    else
        artoffstring=[]; % default: inver arterial aCBV component
    end
else
    artoffstring=[];
end

if nargin >= 9
    if strcmp (spatialoff, 'spatialoff')
        spatialstring= 'off'; % do not perform spatial regularization - much faster for per dynamic CBF analysis for outlierremoval
    else
        spatialstring= 'on'; % default: perform spatial regularization
    end
else
    spatialstring='on';
end

T1t = num2str(SUBJECT.T1t);
T1b = num2str(SUBJECT.T1b);
tau = num2str(SUBJECT.tau);
TR_M0 = num2str(SUBJECT.TR_M0(1));
alpha = num2str(SUBJECT.alpha);
slicetime = num2str(SUBJECT.slicetime/1000); % in s
tic

% JCWS % using BASIL precorrected Mz ->Mxy for LooK Locker; fixbolus, well-mixed fit, , -aif-leadscale=0.00001 as set in BASIL_options.txt (TIP from MARTIN GRAIG, Acknowledge him in a paper!
% note basil option --t1 is the tissue T1 for the quantification, --t1t is the tissue T1 for correction the M0 calibration when TR < 5s
system(['oxford_asl -i ' locationASLlabelcontrolPLDNIFTI ' -c ' locationM0 ' -m ' locationMask ' -o ' outputmap ' ' locationBASILinfostring artoffstring ' --spatial=' spatialstring ' --bolus=' tau ...
    ' --slicedt=' slicetime ' --t1=' T1t ' --t1b=' T1b ' --t1t=' T1t ' --plds=' PLDstring ' --tr=' TR_M0 ' --alpha=' alpha ' --iaf=ct --ibf=tis --casl --fixbolus --cmethod voxel --cgain 1.00']);

disp('BASIL analysis finished')
toc
end
