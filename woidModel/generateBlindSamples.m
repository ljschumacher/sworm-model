% generate random parameter samples
clear all

nSim = 100; % number of samples
nParam = 4; % number of parameters

% seed random number generator for reproducibility
rng('shuffle')

% use rand for uniform prior
samplesRaw = rand(nParam,nSim)'; 
% generate as nParam by nSim and then transpose, so that the random numbers are used
% parameters first, then number of samples. This makes it easier to add
% more samples later on

% scales samples to the appropriate range for the parameters
revRate_range = [0 2];
revRateClusterEdge = samplesRaw(:,1).*(revRate_range(2) - revRate_range(1)) + revRate_range(1);

slowSpeed_range = [0.01 0.33];
slowSpeed = samplesRaw(:,2).*(slowSpeed_range(2) - slowSpeed_range(1)) + slowSpeed_range(1);

Ris_range = [2/3 6];
Ris = samplesRaw(:,3).*(Ris_range(2) - Ris_range(1)) + Ris_range(1);

Rir_range = [2/3 6];
Rir = samplesRaw(:,4).*(Rir_range(2) - Rir_range(1)) + Rir_range(1);

% make a table of the parameters
paramSamples = table(revRateClusterEdge,slowSpeed,Ris,Rir);

% save parameter samples
save(['blindSamples_nSim' num2str(nSim) '_nParam' num2str(nParam)],...
    'paramSamples','revRate_range','slowSpeed_range','Ris_range','Rir_range')