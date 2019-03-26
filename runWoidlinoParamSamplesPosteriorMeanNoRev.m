function [] = runWoidlinoParamSamplesPosteriorMeanNoRev(makeMovie)

% general model parameters for all simulations - unless set otherwise
N = 40; % N: number of objects
M = 18; % M: number of nodes in each object
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
T = 2000;
rc = 0.035; % rc: core repulsion radius (default 0.035 mm)
param.rc = 0;
param.ri = 3*rc;
param.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
param.segmentLength = 1.13/(M - 1);
% -- slow-down parameters --
param.vs = 0.018; % npr1 0.018; N2 0.014
param.slowingNodes = 1:M;% slowingNodes: which nodes register contact (default head and tail)
param.slowingMode = 'stochastic_bynode';
param.k_dwell = 0.0036; % npr1 0.0036; N2 0.25
param.k_undwell = 1.1; % npr1 1.1; N2 0.45
% -- reversal parameters --
param.reversalMode = 'density';
param.revRateClusterEdge = 0;
% -- Lennard-Jones parameters --
param.r_LJcutoff = -1;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
param.sigma_LJ = 0;  % particle size for Lennard-Jones force
param.eps_LJ = 0;
% -- undulation parameters --
param.k_theta = 0;
param.theta_0 = 0;
param.omega_m = 0;
param.deltaPhase = 0;
param.angleNoise = 0.05;
% -- haptotaxis --
param.Rif = 1.2/0.035;
param.haptotaxisMode = 'weighted_additive';
% -- speed and time-step --
param.v0 = [0.33]; % npr1 0.33; N2 0.14
param.dT = min(1/2,rc/param.v0/8); % dT: time step, scales other parameters such as velocities and rates
param.saveEvery = round(1/param.dT);


posteriorfilename = ['../wormtracking/trackingAnalysis/inference/inf_results/'...
    'posteriors_log_PRW_4D_wa_r2_0.01.mat'];
load(posteriorfilename)

% generate samples from posterior
postiSamples = random(posterior{1},1e6);

% truncate boundaries
for dimCtr = 1:size(postiSamples,2)
    overLogIndcs = postiSamples(:,dimCtr)>=supportLimits(2,dimCtr);
    underLogIndcs = postiSamples(:,dimCtr)<=supportLimits(1,dimCtr);
    postiSamples(overLogIndcs|underLogIndcs,:) = [];
end
postiMean = mean(postiSamples);
% set model parameters from posterior
param.drdN_rev = 0;
param.dkdN_dwell = postiMean(2);
param.dkdN_undwell = postiMean(3);
param.f_hapt = 10^postiMean(4);

addpath('visualisation/')

filepath = '~/Dropbox/projects/collectiveBehaviour/sworm-model/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r2/postiPredictiveCheck/';

numReps = 10;
for repCtr = 1:numReps
    filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
        '_v0_' num2str(param.v0) '_vs_' num2str(param.vs) ...
        '_angleNoise_' num2str(param.angleNoise) '_k_theta_' num2str(param.k_theta)...
        '_slow_' param.slowingMode '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell)...
        '_dkdN_' num2str(param.dkdN_dwell,2) '_' num2str(param.dkdN_undwell,2)...
        '_rev' param.reversalMode '_drdN_' num2str(param.drdN_rev,2) ...
        '_haptotaxis_' param.haptotaxisMode '_' num2str(param.f_hapt,2) ...
        '_postiMean_run' num2str(repCtr)];
    
    if ~exist([filepath filename '.mat'],'file')
        rng(repCtr) % for reproducibility
        % run simulation
        disp('Running full-length simulation...')
        [xyarray, currentState, food] = runWoids(T,N,M,L,param);
        xyarray = single(xyarray); % save space by using single precision
        save([filepath filename '.mat'],'xyarray','T','N','M','L','param','currentState')
        
        if makeMovie
            % make movie
            animateWoidTrajectories(xyarray,['movies/woidlinoMovies/paramSampleMovies/' filename '.mp4'],...
                L,0.035);
        end
    end
end