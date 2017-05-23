% run simulations of simplified woid model with single node per woid
% for various speeds, attractions strengths, reversal probabilities...

% issues/todo:
% - always reseed random number generator before each simulation?

clear
close all

% general model parameters for all test - unless set otherwise
N = 204; % N: number of objects
M = 2; % M: number of nodes in each object
L = [20, 20]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
% param.v0 = 1; % v0: speed (default 0.05)
rc = 0.35; % rc: core repulsion radius (default 0.07 mm)
paramAll.rc = 0; % rc: core repulsion radius (default 0.07 mm)
paramAll.segmentLength = 2*rc;
T = 1000; % T: simulation duration
% saveevery = round(1/2/param.dT);
paramAll.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
paramAll.k_l = 1/paramAll.segmentLength; % stiffness of linear springs connecting nodes
% undulations
paramAll.k_theta = 0; % stiffness of rotational springs at nodes
paramAll.omega_m = 0; % angular frequency of oscillation of movement direction, default 0.6 Hz
paramAll.theta_0 = 0; % amplitude of oscillation of movement direction, default pi/4
paramAll.deltaPhase = 0; % for phase shift in undulations and initial positions, default 0.11
% -- reversal parameters --
paramAll.revRate = 0; % revRate: rate for poisson-distributed reversals (default 1/13s)
paramAll.revRateCluster = 0;% revRateCluster: rate for reversals when in a cluster (default 1/130s)
paramAll.revTime = 5; % revTime: duration of reversal events (default 2s, rounded to integer number of time-steps)
paramAll.headNodes = 1;% headNodes: which nodes count as head for defining cluster status, default front 10%
paramAll.tailNodes = 2;% tailNodes: which nodes count as tail for defining cluster status, default back 10%
paramAll.ri = 3*rc;% ri: radius at which worms register contact (default 3 rc)
% -- slow-down parameters --
paramAll.vs = 0;% vs: speed when slowed down (default v0/3)
paramAll.slowingNodes = [];% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
% -- Lennard-Jones parameters --
paramAll.r_LJcutoff = 5*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
% param.eps_LJ = 0;% eps_LJ: strength of LJ-potential
paramAll.sigma_LJ = 2*rc;  % particle size for Lennard-Jones force

paramAll.rc = 0; % turn off contact-forces
revRatesClusterEdge = [0, 0.2, 0.4, 0.6, 0.8];
speeds = [0.15, 0.3];
attractionStrengths = [1e-3, 1e-5, 1e-4];
paramCombis = combvec(revRatesClusterEdge,speeds,attractionStrengths);
nParamCombis = size(paramCombis,2);
parfor paramCtr = 1:nParamCombis
    param = paramAll;
    revRateClusterEdge = paramCombis(1,paramCtr);
    param.revRateClusterEdge = revRateClusterEdge;
    speed = paramCombis(2,paramCtr);
    param.v0 = speed;
    param.dT = min(1/2,rc/param.v0/8); % dT: time step, scales other parameters such as velocities and rates
    saveevery = round(1/2/param.dT);
    attractionStrength = paramCombis(3,paramCtr);
    param.eps_LJ = attractionStrength;
    filename = ['wlM' num2str(M) '_v0_' num2str(param.v0,'%1.0e') '_epsLJ_'...
        num2str(attractionStrength,'%1.0e')...
        '_revRateClusterEdge_' num2str(param.revRateClusterEdge,'%1.0e') '_noContactForces'];
    if ~exist(['results/woidlinos/' filename '.mat'],'file')
        xyarray = runWoids(T,N,M,L,param);
        xyarray = xyarray(:,:,:,1:saveevery:end);
        saveResults(['results/woidlinos/' filename],xyarray,saveevery,T,N,M,L,param)
%         animateWoidTrajectories(xyarray,['tests/woidlinos/' filename],L,rc);
    end
end
