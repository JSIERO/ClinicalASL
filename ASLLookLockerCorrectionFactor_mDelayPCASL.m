function output = ASLLookLockerCorrectionFactor_mDelayPCASL(SUBJECT)
% ClinicalASL toolbox 2023, JCWSiero
% compute Look Locker scaling per PLD
flipangle = SUBJECT.FLIPANGLE;
M0=1;
T1blood= SUBJECT.T1b*1000; % in ms
PLD = SUBJECT.PLDS*1000; %in ms
t = 1:1:PLD(end)*1.1;

deltaPLD=mean(diff(PLD));
T1effblood=1/(1/T1blood-log(cosd(flipangle))/deltaPLD);
Meq_LLblood=M0*(1-exp(-deltaPLD/T1blood))/(1-cosd(flipangle)*exp(-deltaPLD/T1blood)); % Brix et al MRI 1990

Mz_noLL_C = M0*(1 - exp(-t/T1blood)); % control data
Mz_noLL_L = M0*(1 - 2* exp(-t/T1blood)); % label data
deltaM_noLL = M0*exp(-t/T1blood);

t1=1:PLD(1);
Mz_LL_C(t1) = M0*(1 - exp(-t1/T1blood));
Mz_LL_L(t1) = M0*(1 - 2* exp(-t1/T1blood));
t2=PLD(1)+1:t(end);
Mz_LL_C(t2) = Mz_LL_C(PLD(1))*exp(-t(t2-PLD(1))/T1effblood) + Meq_LLblood*(1 - exp(-t(t2-PLD(1))/T1effblood));
Mz_LL_L(t2) = Mz_LL_L(PLD(1))*exp(-t(t2-PLD(1))/T1effblood) + Meq_LLblood*(1 - exp(-t(t2-PLD(1))/T1effblood));
deltaM_LL = Mz_LL_C - Mz_LL_L;

deltaM_ratio_LL_noLL=deltaM_LL./deltaM_noLL;

tpointsMxy = PLD;

deltaMxy_ratio_LL_noLL=  deltaM_ratio_LL_noLL*sind(flipangle);
deltaMxy_ratio_LL_noLL_forBASIL = deltaMxy_ratio_LL_noLL(tpointsMxy);
disp([' Look-Locker correction factor for flipangle = ' num2str(flipangle) '(deg), PLDs(ms) :' num2str(PLD) ', and deltaPLD(ms) = ' num2str(deltaPLD) ':  ' num2str(deltaMxy_ratio_LL_noLL_forBASIL)])

output = deltaMxy_ratio_LL_noLL_forBASIL;