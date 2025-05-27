 function ASLQASL_SSVBanalysis(SUBJECT, locationASLlabelcontrolPLDNIFTI, locationM0, locationMask, outputmap, PLDlist, artoff, spatialoff)
% ClinicalASL toolbox 2023, JCWSiero
PLDstring = sprintf('%.05g,' , SUBJECT.PLDS(PLDlist));
PLDstring = PLDstring(1:end-1);% strip final comma

if nargin >= 7
    if strcmp (artoff, 'artoff')
        artoffstring= ' --artoff'; % do not inver arterial aCBV component, needed for aterial transit time aterfact (ATA) mapping - also much faster, use for per dynamic CBF analysis (outlierremoval)
    else
        artoffstring=[]; % default: inver arterial aCBV component
    end
else
    artoffstring=[];
end

if nargin >= 8
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

% Build the Docker run command for SSVB method
cmd = sprintf([...
    'docker run -u $(id -u):$(id -g) --rm -v "%s":/input -v "%s":/output ', ...
    'quantifiedimaging/qasl ', ...
    '-i /input/' locationASLlabelcontrolPLDNIFTI ' ', ...
    '-o /output/' outputmap ' ', ...
    '-c /input/' locationM0 ' ', ...
    '-m /input/' locationMask ' ', ...
    artoffstring ' ',...
    spatialstring ' ',...
    '--inference-method=ssvb ', ...
    '--tau=' tau ' --slicedt=' slicetime ' --t1=' T1t ' --t1b=' T1b ' --t1t=' T1t ' ', ...
    '--plds=' PLDstring ' ', ...
    '--tr=' TR_M0 ' --alpha=' alpha ' --iaf=ct --ibf=tis --casl --cgain=1.0 --readout=2D --overwrite --save-calib'], ...
    SUBJECT.ASLdir, SUBJECT.ASLdir);

% Run the Dockerized QASL
status = system(cmd);
if status ~= 0
    error('Quantified Imaging QASL docker command failed.');
end

toc
end
