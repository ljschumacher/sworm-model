% generate random parameter samples
clear all

nSamples = 1e5; % number of samples
nParam = 4; % number of parameters

% seed random number generator for reproducibility
rng(2)

% load posterior gmmodel
posteriorfilename = ['../../wormtracking/trackingAnalysis/inference/inf_results/'...
    'posteriors_log_PRW_4D_wa_r1_0.1.mat'];
load(posteriorfilename)

posti_N2 = posterior{2};
% sample from prior
samplesRaw = random(posti_N2,nSamples);
% adjust support limits
supportLimits(1,4) = -4;
supportLimits(:,4) = 10.^supportLimits(:,4); % adjust for log10-basis parameters
% enforce support Limits
while any(any(samplesRaw<=supportLimits(1,:) |...
        samplesRaw>=supportLimits(2,:),2))
    replaceLogIdcs = any(samplesRaw<=supportLimits(1,:) |...
        samplesRaw>=supportLimits(2,:),2);
    nReplace = nnz(replaceLogIdcs);
    samplesRaw(replaceLogIdcs,:) = random(posti_N2,nReplace);
end

drdN_rev = samplesRaw(:,1);

dkdN_dwell = samplesRaw(:,2);

dkdN_undwell = samplesRaw(:,3);

f_hapt = samplesRaw(:,4);

% make a table of the parameters
paramSamples = table(drdN_rev,dkdN_dwell,dkdN_undwell,f_hapt);

% save parameter samples
save(['paramSamples_nSamples' num2str(nSamples) '_log_PRW_4D_wa_r2_N2' ...
    '.mat'],'paramSamples','supportLimits')