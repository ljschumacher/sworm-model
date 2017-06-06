% run simulations of simplified woid model with single node per woid
% for various speeds, attractions strengths, reversal probabilities...

% issues/todo:
% - always reseed random number generator before each simulation?

clear
close all

% general model parameters for all test - unless set otherwise
N = 40; % N: number of objects
M = 49; % M: number of nodes in each object
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
T = 500; % T: simulation duration
rc = 0.035;
% saveevery = round(1/2/param.dT);
paramAll.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
% -- slow-down parameters --
paramAll.vs = 0;% vs: speed when slowed down (default v0/3)
paramAll.slowingNodes = [1:max(round(M/10),1) (M-max(round(M/10),1)+1):M];% slowingNodes: which nodes register contact (default head and tail)
% -- Lennard-Jones parameters --
paramAll.r_LJcutoff = 5*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
paramAll.sigma_LJ = 2*rc;  % particle size for Lennard-Jones force

revRatesClusterEdge = [0, 0.1, 0.2, 0.4, 0.8];
speeds = [0.33];
slowspeeds = [0.33, 0.1, 0.05, 0.025];
attractionStrengths = [0];
paramCombis = combvec(revRatesClusterEdge,speeds,slowspeeds,attractionStrengths);
nParamCombis = size(paramCombis,2);
parfor paramCtr = 1:nParamCombis
    param = paramAll;
    revRateClusterEdge = paramCombis(1,paramCtr);
    param.revRateClusterEdge = revRateClusterEdge;
    speed = paramCombis(2,paramCtr);
    param.omega_m = 2*pi*0.6/0.33*speed;
    param.v0 = speed;
    param.dT = min(1/2,rc/param.v0/16); % dT: time step, scales other parameters such as velocities and rates
    saveevery = round(1/4/param.dT);
    param.vs = paramCombis(3,paramCtr);
    attractionStrength = paramCombis(4,paramCtr);
    if attractionStrength>0
        param.r_LJcutoff = 5*rc;
    else
        param.r_LJcutoff = -1; % don't need to compute attraction if it's zero
    end
    param.eps_LJ = attractionStrength;
    filename = ['woids_v0_' num2str(param.v0,'%1.0e') ...
        '_vs_' num2str(param.vs,'%1.0e') ...
        '_gradualSlowDown_' ...
        '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
        '_revRateClusterEdge_' num2str(param.revRateClusterEdge,'%1.0e')];
    if ~exist(['results/woids/' filename '.mat'],'file')
        rng(1) % set random seed to be the same for each simulation
        xyarray = runWoids(T,N,M,L,param);
        xyarray = xyarray(:,:,:,1:saveevery:end);
        saveResults(['results/woids/' filename],...
        struct('xyarray',xyarray,'saveevery',saveevery,'T',T,'N',N,'M',M,'L',L,'param',param))
    end
end
