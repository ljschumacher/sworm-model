% run simulations of simplified woid model with single node per woid
% for various speeds, attractions strengths, reversal probabilities...

% issues/todo:
% - always reseed random number generator before each simulation?

clear
close all

% general model parameters for all test - unless set otherwise
N = 204; % N: number of objects
M = 1; % M: number of nodes in each object
L = [20, 20]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
% param.v0 = 1; % v0: speed (default 0.05)
rc = 0.5; % rc: core repulsion radius (default 0.07 mm)
% param.segmentLength = 0;
% param.dT = param.rc/param.v0/4; % dT: time step, scales other parameters such as velocities and rates
% T = 1000/param.dT; % T: simulation duration (number of time-steps)
% saveevery = round(1/2/param.dT);
param.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
param.k_l = 0; % stiffness of linear springs connecting nodes
% undulations
param.k_theta = 0; % stiffness of rotational springs at nodes
param.omega_m = 0; % angular frequency of oscillation of movement direction, default 0.6 Hz
param.theta_0 = 0; % amplitude of oscillation of movement direction, default pi/4
param.deltaPhase = 0; % for phase shift in undulations and initial positions, default 0.11
% -- reversal parameters --
param.revRate = 0.1; % revRate: rate for poisson-distributed reversals (default 1/13s)
param.revRateCluster = 0.1;% revRateCluster: rate for reversals when in a cluster (default 1/130s)
param.revTime = 2; % revTime: duration of reversal events (default 2s, rounded to integer number of time-steps)
param.headNodes = 1;% headNodes: which nodes count as head for defining cluster status, default front 10%
param.tailNodes = [];% tailNodes: which nodes count as tail for defining cluster status, default back 10%
param.ri = 3*rc;% ri: radius at which worms register contact (default 3 rc)
% -- slow-down parameters --
param.vs = 0;% vs: speed when slowed down (default v0/3)
param.slowingNodes = [];% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
% -- Lennard-Jones parameters --
param.r_LJcutoff = 5*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
% param.eps_LJ = 0;% eps_LJ: strength of LJ-potential
param.sigma_LJ = 2*rc;  % particle size for Lennard-Jones force

param.rc = 0; % turn off contact-forces
for revRate = [0, 0.1, 1]
    param.revRate = revRate;
    param.revRateCluster = revRate;
    for speed = [0.1, 0.5, 1]
        param.v0 = speed;
        param.dT = min(1/2,rc/param.v0/4); % dT: time step, scales other parameters such as velocities and rates
        T = 1000; % T: simulation duration
        saveevery = round(1/2/param.dT);
        for attractionStrength = [1e-5, 1e-4, 5e-4]
            param.eps_LJ = attractionStrength;
            filename = ['wl_v0_' num2str(param.v0,'%1.0e') '_epsLJ_' num2str(attractionStrength,'%1.0e')...
                '_revRate_' num2str(param.revRate,'%1.0e') '_noContactForces'];
            if ~exist(['results/woidlets/' filename '.mat'],'file')
                xyarray = runWoids(T,N,M,L,param);
                xyarray = xyarray(:,:,:,1:saveevery:end);
                save(['results/woidlets/' filename])
                animateWoidTrajectories(xyarray,['tests/woidlets/' filename],L,rc);
            end
        end
    end
end