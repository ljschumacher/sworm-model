% generate random parameter samples

nSim = 1000; % number of samples
nParam = 4; % number of parameters

% seed random number generator for reproducibility
rng(1)

% use rand for uniform prior
samplesRaw = rand(nParam,nSim)'; 
% generate as nParam by nSim and then transpose, so that the random numbers are used
% parameters first, then number of samples. This makes it easier to add
% more samples later on

% scales samples to the appropriate range for the parameters
revRateRange = [0 2];
revRateClusterEdge = samplesRaw(:,1).*(revRateRange(2) - revRateRange(1)) + revRateRange(1);

slowSpeedRange = [0.01 0.33];
slowSpeed = samplesRaw(:,2).*(slowSpeedRange(2) - slowSpeedRange(1)) + slowSpeedRange(1);

Ris_Range = [2 6];
Ris = samplesRaw(:,3).*(Ris_Range(2) - Ris_Range(1)) + Ris_Range(1);

Rir_Range = [2 6];
Rir = samplesRaw(:,3).*(Rir_Range(2) - Rir_Range(1)) + Rir_Range(1);

% make a table of the parameters
paramSamples = table(revRateClusterEdge,slowSpeed,Ris,Rir);

% save parameter samples
save(['paramSamples_nSim' num2str(nSim) '_nParam' num2str(nParam)])