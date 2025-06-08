function T1fromPLD=ASLT1fromPLDCompute(DATA4D, MASK, TIMEARRAY)
% ClinicalASL toolbox 2023, JCWSiero
% compute T1w image from multi-delay PCASL (Look-Locker) PLD (not the PLD) data across the different PLDs
% we observe that the PLD across PLDs show a exponential decay that depends on T1
% take natural log and fit for 'T1'

dims = size(DATA4D);
DATA2D = reshape(DATA4D, dims(1)*dims(2)*dims(3),dims(4));

%take the natural log of data
%DATA2D = log(DATA2D); 

CONSTANTARRAY = ones(size(DATA2D, 2), 1); % create array with ones, for the constant term
T1fit = zeros(size(DATA2D,1),2);
warning off
for i = 1:size(DATA2D,1)
A = [CONSTANTARRAY, DATA2D(i,:)']; % construct system of linear equations: Ax=b
T1fit(i,:) = A\TIMEARRAY'; % solve the system of linear equations, yielding the decay rate constant (R1 weighted) and offset (not used)
end
warning on

data_R1fit = reshape(T1fit(:,2),dims(1),dims(2),dims(3)); % reshape to construct a 3D map of the decay rate constant from the linear fit
data_T1fit_brain = data_R1fit.*MASK*1e3; % convert decay rate constant to the decay time constant (T1 weighted) and impose mask to include brain-only voxels
data_T1fit_brain(isnan(data_T1fit_brain)) = 0;

data_T1fit_brain_mask_high = data_T1fit_brain<=2; % mask to exclude too high values after fit 
data_T1fit_brain_mask_low = data_T1fit_brain>0; % mask to exclude negative values after fit (e.g. CSF)

T1fromPLD = data_T1fit_brain.*data_T1fit_brain_mask_high.*data_T1fit_brain_mask_low;
T1fromPLD(isnan(T1fromPLD)) = 0;
