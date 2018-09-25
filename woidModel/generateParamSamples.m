% generate random parameter samples
clear all

nSamples = 1e5; % number of samples
nParam = 4; % number of parameters

% seed random number generator for reproducibility
rng(1)

% set parameters
M = 18;
angleNoise = 0.05;
k_theta = 0;
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
reversalMode = 'density';
haptotaxisMode = 'weighted_additive';

% load reduced prior gmmodel
load(['priors4D_M_' num2str(M) '_noVolExcl' ...
    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
    '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
    '_rev' reversalMode '_haptotaxis_' haptotaxisMode ...
    '.mat'],'supportLimits','prior_npr1')

% sample from prior
samplesRaw = random(prior_npr1,nSamples);

% enforce support Limits
while any(any(samplesRaw<=supportLimits(1,:) |...
        samplesRaw>=supportLimits(2,:),2))
    replaceLogIdcs = any(samplesRaw<=supportLimits(1,:) |...
        samplesRaw>=supportLimits(2,:),2);
    nReplace = nnz(replaceLogIdcs);
    samplesRaw(replaceLogIdcs,:) = random(prior_npr1,nReplace);
end

drdN_rev = samplesRaw(:,1);

dkdN_dwell = samplesRaw(:,2);

dkdN_undwell = samplesRaw(:,3);

f_hapt = samplesRaw(:,4);

% make a table of the parameters
paramSamples = table(drdN_rev,dkdN_dwell,dkdN_undwell,f_hapt);

% save parameter samples
save(['paramSamples_nSamples' num2str(nSamples) '_nParam' num2str(nParam) ...
    '_M_' num2str(M) '_noVolExcl' ...
    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
    '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
    '_rev' reversalMode '_haptotaxis_' haptotaxisMode ...
    '.mat'],'paramSamples','supportLimits')