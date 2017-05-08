% test simplified woid model with single node per woid

% issues/todo:
% - always reseed random number generator before each simulation?

clear
close all

% general model parameters for all test - unless set otherwise
N = 160; % N: number of objects
M = 1; % M: number of nodes in each object
L = 10;%[20, 20]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
param.v0 = 1; % v0: speed (default 0.05)
param.rc = 0.5; % rc: core repulsion radius (default 0.07 mm)
% param.segmentLength = 0;
param.dT = param.rc/param.v0/4; % dT: time step, gets adapted in simulation
T = 100; % T: simulation duration (number of time-steps)
saveevery = round(1/2/param.dT);
param.bc = 'noflux'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
param.k_l = 0; % stiffness of linear springs connecting nodes
% undulations
param.k_theta = 0; % stiffness of rotational springs at nodes
param.omega_m = 0; % angular frequency of oscillation of movement direction, default 0.6 Hz
param.theta_0 = 0; % amplitude of oscillation of movement direction, default pi/4
param.deltaPhase = 0; % for phase shift in undulations and initial positions, default 0.11
% -- reversal parameters --
param.revRate = 0; % revRate: rate for poisson-distributed reversals (default 1/13s)
param.revRateCluster = 0;% revRateCluster: rate for reversals when in a cluster (default 1/130s)
param.revTime = 0; % revTime: duration of reversal events (default 2s, rounded to integer number of time-steps)
param.headNodes = 1;% headNodes: which nodes count as head for defining cluster status, default front 10%
param.tailNodes = [];% tailNodes: which nodes count as tail for defining cluster status, default back 10%
param.ri = 3*param.rc;% ri: radius at which worms register contact (default 3 rc)
% -- slow-down parameters --
param.vs = param.v0;% vs: speed when slowed down (default v0/3)
param.slowingNodes = [];% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
% -- Lennard-Jones parameters --
param.r_LJcutoff = param.ri;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
param.eps_LJ = 0;% eps_LJ: strength of LJ-potential

xyarray = runWoids(T,N,M,L,param);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/woidlets/test_noflux_circular',L,param.rc);

L = [20, 20];
N = 204;
xyarray = runWoids(T,N,M,L,param);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/woidlets/test_noflux_square',L,param.rc);

param.bc = 'periodic';
xyarray = runWoids(T,N,M,L,param);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/woidlets/test_periodic_square',L,param.rc);