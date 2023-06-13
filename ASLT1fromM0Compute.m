function T1fromM0=ASLT1fromM0Compute(DATA4D, MASK, TIMEARRAY)
% ClinicalASL toolbox 2023, JCWSiero
% compute T1w image from multi-delay PCASL (Look-Locker) M0 data across the different PLDs
% we observe the M0 accros PLDs show a exponential decay  that depends on T1
% take natural log and fit for 'T1'

dims = size(DATA4D);
DATA2D = reshape(DATA4D, dims(1)*dims(2)*dims(3),dims(4));

%take natural log of data
DATA2D = log(DATA2D);

CONSTANTARRAY = ones(size(DATA2D, 2), 1); % create array with ones, for the constant term
T1fit = zeros(size(DATA2D,1),2);
warning off
for i = 1:size(DATA2D,1)
A = [CONSTANTARRAY, DATA2D(i,:)'];
T1fit(i,:) = A\TIMEARRAY';
end
warning on

data_T1fit2 = reshape(T1fit(:,2),dims(1),dims(2),dims(3));
data_T1fit2_brain = -1/data_T1fit2.*MASK*1e3;
data_T1fit2_brain(isnan(data_T1fit2_brain)) = 0;
data_T1fit2_brain_mask_high = data_T1fit2_brain<=300; % mask to exlude too high values after fit 
data_T1fit2_brain_mask_low = data_T1fit2_brain>0; % mask to exlude negative values after fit (CSF)

T1fromM0 = data_T1fit2_brain.*data_T1fit2_brain_mask_high.*data_T1fit2_brain_mask_low;
T1fromM0(isnan(T1fromM0)) = 0;