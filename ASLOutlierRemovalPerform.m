function [Ioutlier_allsteps, Ioutlier_step1, Ioutlier_step2, NoOutliers_logical] = ASLOutlierRemovalPerform(CBFdata4D, brainmask, outlierFactor, GMmask, WMmask, CSFmask)
%% ASL OutlierRemoval (volumes), Siero 2019 - 01-03-2019
% SCORE method by Duloi et al JMRI 2017
% now works with NaN

%confine to brain tissue voxels
GMmask = logical(GMmask.*brainmask);
WMmask = logical(WMmask.*brainmask);
CSFmask = logical(CSFmask.*brainmask);
mean_cbf_gm = [];

for i = 1:size(CBFdata4D,4)
  dummy = squeeze(CBFdata4D(:,:,:,i));
  mean_cbf_gm(i) = mean(dummy(GMmask),'omitnan');
  if isnan(mean_cbf_gm(i))
    error(['NaN value found in meanCBF volume : ' num2str(i)])
  elseif isinf(mean_cbf_gm(i))
    disp(['Inf value found in meanCBF volume : ' num2str(i)])
    disp(['Setting Inf values to zero for ' num2str(nnz(isinf(dummy))) ' voxels']);
    dummy(isinf(dummy)) = 0;
  end
  mean_cbf_gm(i) = mean(dummy(GMmask),'omitnan');
  CBFdata4D(:,:,:,i) = dummy;
end

disp('************************************************** Perform Ooutlier Removal SCORE by Duloi et al **********************************************')
Ioutlier_step1=[];
Ioutlier_step2=[];
mu_cbf = median(mean_cbf_gm);
std_cbf = 1.4826*mad(mean_cbf_gm,1); %median absolute deviation

step1 = abs(mean_cbf_gm - mu_cbf)<= std_cbf*outlierFactor; % creates logical array of length NREPEATS, with 0's where CBF > threshold
data_step1 = CBFdata4D(:,:,:,step1);
if find(step1 == 0)
    disp(['Step1: Volume(s) removed (|meanCBF - mu| > ' num2str(outlierFactor) '*std): ' num2str(find(step1 == 0))]);    
else
    disp(['Step1: Volume(s) removed (|meanCBF - mu| > ' num2str(outlierFactor) '*std): 0']);
end
disp('')
Ioutlier_step1 = find(step1 == 0);
Ioutlier_allsteps = Ioutlier_step1; % initialize Ioutlier_allsteps
mean_cbf_step2 = mean(data_step1,4,'omitnan');

var_gm = var(mean_cbf_step2(GMmask));
N_gm = length(mean_cbf_step2(GMmask))-1;%unbiased version for pooled variance

if nargin > 4
  var_csf = var(mean_cbf_step2(CSFmask),'omitnan');
  var_wm = var(mean_cbf_step2(WMmask),'omitnan');
  N_csf = length(mean_cbf_step2(CSFmask))-1;
  N_wm = length(mean_cbf_step2(WMmask))-1;
  V = (N_gm*var_gm + N_csf*var_csf + N_wm*var_wm)/(N_gm + N_csf + N_wm); %pooled variance
  disp(' Using Gray Matter, White Matter and CSF Spatial Pooled Variance')
else
  V = var_gm; % variance
  disp(' Using Gray Matter, Spatial Variance')
end

Vprev = V + 0.01*V;
data_step2 = data_step1;

mean_cbf_prev_mask = mean_cbf_step2(brainmask);
disp(['Spatial Pooled Variance, iter = 0 : ' num2str(V)]);
iter = 0;
while V < Vprev && size(data_step2,4) > 1
  Vprev = V;
  cbfvols_inmask = [];
  iter = iter+1;
  for i = 1:size(data_step2,4)
    dummy = data_step2(:,:,:,i);
    cbfvols_inmask(:,i) = dummy(brainmask);
  end
  R = zeros(1,size(cbfvols_inmask,2));
  for i = 1:size(cbfvols_inmask,2)
    Rarray = corrcoef(cbfvols_inmask(:,i),mean_cbf_prev_mask, 'rows','complete'); % works with NaN entries
    R(i) = Rarray(2);
  end
  data_step3 = data_step2;
  dummy2 = data_step3(:,:,:,R == max(R));
  meancbf_removed = mean(dummy2(GMmask),'omitnan');
  data_step3(:,:,:,R == max(R)) = [];% throw out volume with max corr, but only if spatial variance increases
  mean_cbf_prev = mean(data_step3,4,'omitnan');
  mean_cbf_prev_mask = mean_cbf_prev(brainmask);
  var_gm = var(mean_cbf_prev(GMmask),'omitnan');
  N_gm = length(mean_cbf_prev(GMmask))-1;%unbiased version for pooled variance

  if nargin > 4
    var_csf = var(mean_cbf_prev(CSFmask),'omitnan');
    var_wm = var(mean_cbf_prev(WMmask),'omitnan');
    N_csf = length(mean_cbf_prev(CSFmask))-1;
    N_wm = length(mean_cbf_prev(WMmask))-1;
    V = (N_gm*var_gm + N_csf*var_csf + N_wm*var_wm)/(N_gm + N_csf + N_wm); %pooled variance
  else
    V = var_gm; % variance
  end
  if V < Vprev
    dummy2 = data_step2(:,:,:,R == max(R));
    meancbf_removed = mean(dummy2(GMmask),'omitnan');
    Ioutlier_step2 = find(meancbf_removed == mean_cbf_gm);
    data_step2 = data_step3;
    disp(['Spatial Pooled Variance, iter = ' num2str(iter) ' : ' num2str(V)]);
    disp(['Step2: Volume(i) removed = ' num2str(Ioutlier_step2)]);
    Ioutlier_allsteps = sort([Ioutlier_allsteps Ioutlier_step2]);
  else
    disp(['Spatial Pooled Variance, iter = ' num2str(iter) ' : ' num2str(V)]);
    disp(['Step2: Pooled Variance Increased, STOPPED ']);
  end
end

NoOutliers_logical(1:size(CBFdata4D,4)) = 1;
NoOutliers_logical(Ioutlier_allsteps) = 0;
NoOutliers_logical = logical(NoOutliers_logical);
disp(['Outlier volumes all steps (total number): ' num2str(Ioutlier_allsteps) ' (' num2str(numel(Ioutlier_allsteps)) ')'])
disp('*************************************************** CBF Outlier identification finished ***********************************');
