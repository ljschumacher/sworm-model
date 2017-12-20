function [] = runWoidlinoParamSamples(sampleCtr)
% run simulations of simplified woid model with single node per woid
% for previously generated random parameter samples

% general model parameters for all simulations - unless set otherwise
N = 40; % N: number of objects
M = 18; % M: number of nodes in each object
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
T = 1000;
rc = 0.035; % rc: core repulsion radius (default 0.035 mm)
param.rc = rc;
param.ri = 3*rc;
param.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
param.segmentLength = 1.13/(M - 1);
param.k_l = 40; % stiffness of linear springs connecting nodes
% -- slow-down parameters --
param.vs = 0.018;
param.slowingNodes = 1:M;% slowingNodes: which nodes register contact (default head and tail)
param.slowingMode = 'stochastic_bynode';
param.k_dwell = 0.0036; 
param.k_undwell = 1.1;
% -- Lennard-Jones parameters --
param.r_LJcutoff = 4*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
param.sigma_LJ = 2*rc;  % particle size for Lennard-Jones force
param.eps_LJ = 0;
if param.eps_LJ>0
    param.r_LJcutoff = 5*rc;
else
    param.r_LJcutoff = -1; % don't need to compute attraction if it's zero
end
% -- undulation parameters --
param.theta_0 = 0;
param.omega_m = 0;
param.deltaPhase = 0;
% -- speed and time-step --
param.v0 = [0.33];
param.dT = min(1/2,rc/param.v0/16); % dT: time step, scales other parameters such as velocities and rates
param.saveEvery = round(1/2/param.dT);

% load randomly generated parameter samples
load('paramSamples_nSim1000_nParam2.mat','paramSamples')
% set model parameters from generated samples
param.revRateClusterEdge = paramSamples.revRateClusterEdge(sampleCtr);
param.dkdN_dwell = paramSamples.dkdN(sampleCtr);

filename = ['/work/lschumac/woidlinos/wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
    '_v0_' num2str(param.v0) '_vs_' num2str(param.vs) ...
    ...'_Ris_' num2str(param.Ris) ...
    '_' param.slowingMode 'SlowDown' '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
    '_dkdN_' num2str(param.dkdN_dwell)...
    '_revRateClusterEdge_' num2str(param.revRateClusterEdge) ...
    ...'_Rir_' num2str(param.Rir)
    '_sample_' num2str(sampleCtr)];
if ~exist([filename '.mat'],'file')
    rng('shuffle') % set random seed to be DIFFERENT for each simulation
    [xyarray, currentState] = runWoids(T,N,M,L,param);
    xyarray = single(xyarray); % save space by using single precision
    save([filename '.mat'],'xyarray','T','N','M','L','param','currentState')
end
end