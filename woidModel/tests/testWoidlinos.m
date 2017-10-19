% test simplified woid model with two nodes per woid

% issues/todo:
% - always reseed random number generator before each simulation?

clear
close all

addpath('../')
addpath('../visualisation')
% general model parameters for all test - unless set otherwise
M = 2; % M: number of nodes in each object
L = 10;%[20, 20]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
param.v0 = 0.5; % v0: speed (default 0.05)
rc = 0.35;
param.rc = -1; % rc: core repulsion radius (default 0.07 mm)
param.segmentLength = 2*rc;
param.dT = rc/param.v0/8; % dT: time step, gets adapted in simulation
T = 50; % T: simulation duration (number of time-steps)
param.saveEvery = round(1/2/param.dT);
param.k_l = 40; % stiffness of linear springs connecting nodes
param.bc = 'noflux'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
% undulations
param.k_theta = 0; % stiffness of rotational springs at nodes
param.omega_m = 0; % angular frequency of oscillation of movement direction, default 0.6 Hz
param.theta_0 = 0; % amplitude of oscillation of movement direction, default pi/4
param.deltaPhase = 0; % for phase shift in undulations and initial positions, default 0.11
% -- reversal parameters --
% param.revRate = 0; % revRate: rate for poisson-distributed reversals (default 1/13s)
% param.revRateCluster = 0;% revRateCluster: rate for reversals when in a cluster (default 1/130s)
% param.revTime = 2; % revTime: duration of reversal events (default 2s, rounded to integer number of time-steps)
% param.headNodes = [];% headNodes: which nodes count as head for defining cluster status, default front 10%
% param.tailNodes = [];% tailNodes: which nodes count as tail for defining cluster status, default back 10%
param.ri = 3*rc;% ri: radius at which worms register contact (default 3 rc)
% -- slow-down parameters --
param.vs = param.v0;% vs: speed when slowed down (default v0/3)
param.slowingNodes = [];% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
% -- Lennard-Jones parameters --
param.r_LJcutoff = param.ri;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
param.eps_LJ = 0;% eps_LJ: strength of LJ-potential

% N = 20;
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_noflux_circular',L,rc);

% L = [20, 20];
% N = 50;
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_noflux_square',L,rc);

param.bc = 'periodic';
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_periodic_square',L,rc);

% param.r_LJcutoff = 2*rc;
% param.eps_LJ = 1e-2;
% param.sigma_LJ = 2*rc;
% rng(1)
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_periodic_square_repulsionOnly',L,rc);
% param.eps_LJ = 0;
% 
% param.r_LJcutoff = rc;
% param.rc = rc;
% rng(1)
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_periodic_square_contactForces',L,rc);
% 
% param.revRate = 0.1;
% param.revTime = 10;
% N = 1;
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_periodic_square_reversals',L,rc);
% 
% param.revRate = 0;
% param.revRateClusterEdge = 1;
% param.revTime = 5;
% param.headNodes = 1;
% param.tailNodes = 2;
% N = 50;
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_periodic_square_reversalsClusterEdge',L,rc);

% % test woidlinos with adhesion and wout reversals/slowing
% rng(1)
% T= 400;
% N = 40;
% L = [7.5, 7.5];
% param.rc = 0.105;
% param.ri = 3*param.rc;
% param.segmentLength = 2*param.rc;
% param.v0 = 0.33;
% param.dT = min(1/2,param.rc/param.v0/8);
% param.saveEvery = round(1/4/param.dT);
% param.vs = param.v0;
% param.revRate = 0;
% param.revRateClusterEdge = 0;
% param.revRateCluster = 0;
% param.r_LJcutoff = 5*param.rc;
% param.eps_LJ = 1e-3;
% param.sigma_LJ = 2*param.rc;
% param.LJnodes = 1;
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,...
%     ['woidlino_test_movies/test_periodic_square_adhesion_eps'...
%     num2str(param.eps_LJ,'%1.0e') '_noRev'...
%     '_LJhead'],L,param.rc);

% % test long woidlinos without volume exclusion
rng(1)
T= 40;
N = 40;
M = 18;
L = [7.5, 7.5];
rc0 = 0.035;
param.bc = 'periodic';
param.rc = 0;
param.segmentLength = (1.2-2*rc0)/(M - 1);
param.v0 = 0.33;
param.dT = min(1/2,rc0/param.v0/8);
param.saveEvery = round(1/4/param.dT);
param.vs = param.v0;
param.revRate = 0;
param.revRateClusterEdge = 0;
param.revRateCluster = 0;
param.r_LJcutoff = 5*rc0;
param.eps_LJ = 0;
param.sigma_LJ = 2*rc0;
param.LJnodes = 1:M;
param.angleNoise = 0.01;
xyarray = runWoids(T,N,M,L,param);
animateWoidTrajectories(xyarray,...
    ['woidlino_test_movies/test_longBody_periodic_square_noRev'...
    'eps_LJ_' num2str(param.eps_LJ,'%1.0e'),...
    '_noVolExcl' '_angleNoise'],L,rc0);