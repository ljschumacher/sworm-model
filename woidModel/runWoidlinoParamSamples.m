function [] = runWoidlinoParamSamples(sampleCtr)
% run simulations of simplified woid model with single node per woid
% for previously generated random parameter samples

% general model parameters for all simulations - unless set otherwise
N = 40; % N: number of objects
M = 2; % M: number of nodes in each object
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
T = 500;
rc = 0.105; % rc: core repulsion radius (default 0.035 mm)
param.rc = rc;
param.ri = 3*rc;
% saveevery = round(1/2/param.dT);
param.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
param.segmentLength = 2*rc;%1.2/(M - 1);
param.k_l = 40; % stiffness of linear springs connecting nodes
% -- slow-down parameters --
param.slowingNodes = 1:M;% slowingNodes: which nodes register contact (default head and tail)
param.headNodes = 1:max(round(M/10),1);
param.tailNodes = (M-max(round(M/10),1)+1):M;
% -- Lennard-Jones parameters --
param.r_LJcutoff = 5*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
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
param.saveEvery = round(1/4/param.dT);

% load randomly generated parameter samples
load('paramSamples_nSim1000_nParam4.mat','paramSamples')
% set model parameters from generated samples
param.revRateClusterEdge = paramSamples.revRateClusterEdge(sampleCtr);
param.vs = paramSamples.slowSpeed(sampleCtr);
param.Rir = paramSamples.Rir(sampleCtr);
param.Ris = paramSamples.Ris(sampleCtr);

filename = ['/work/lschumac/woidlinos/wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
    '_v0_' num2str(param.v0) '_vs_' num2str(param.vs) ...
    '_Ris_' num2str(param.Ris) ...
    '_revRateClusterEdge_' num2str(param.revRateClusterEdge) ...
    '_Rir_' num2str(param.Rir) '_sample_' sampleCtr];
if ~exist([filename '.mat'],'file')
    rng(1) % set random seed to be the same for each simulation
    xyarray = runWoids(T,N,M,L,param);
    xyarray = single(xyarray); % save space by using single precision
    save([filename '.mat'],'xyarray','T','N','M','L','param')
end
end