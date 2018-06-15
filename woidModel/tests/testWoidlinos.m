% test simplified woid model with rods

% issues/todo:
% - always reseed random number generator before each simulation?

clear
close all

addpath('../')
addpath('../visualisation')
addpath('../analysis')

% general model parameters for all test - unless set otherwise
M = 18; % M: number of nodes in each object
L = [7.5, 7.5];%[20, 20]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
param.v0 = 0.33; % v0: speed (default 0.05)
rc = 0.035;
param.rc = 0; % rc: core repulsion radius (default 0.07 mm)
param.segmentLength = 1.13/(M - 1);
param.dT = rc/param.v0/16; % dT: time step, gets adapted in simulation
T = 50; % T: simulation duration (number of time-steps)
param.saveEvery = round(1/2/param.dT);
param.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
% undulations
param.omega_m = 0; % angular frequency of oscillation of movement direction, default 0.6 Hz
param.theta_0 = 0; % amplitude of oscillation of movement direction, default pi/4
param.deltaPhase = 0; % for phase shift in undulations and initial positions, default 0.11
% -- reversal parameters --
param.ri = 3*rc;% ri: radius at which rods register contact (default 3 rc)
% -- slow-down parameters --
param.vs = param.v0;% vs: speed when slowed down (default v0/3)
param.slowingNodes = [];% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
% -- Lennard-Jones parameters --
param.r_LJcutoff = 4*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
param.sigma_LJ = 2*rc;
param.eps_LJ = 0;% eps_LJ: strength of LJ-potential
food = [];

% % test clustered initial conditions
% param.bc = 'free';
% param.sigma_LJ = 0;
% L = [3.6, 3.6];
% xyarray = runWoids(5,40,18,L,param);
% animateWoidTrajectories(xyarray,['woidlino_test_movies/test_clustered'],L);

% % test angle noise 
% param.angleNoise = 0.02;% not much point making this any bigger than 10, because it's angular
% param.bc = 'free';
% param.k_theta = 0;
% L = [3 3];
% rng(1)
% xyarray = runWoids(100,1,M,L,param);
% filename = ['woidlino_test_movies/test_free_'...
%     'angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)];

% % test haptotaxis
% rng(2)
% L = [3 3];
% param.f_hapt = 0.25;
% xyarray = runWoids(20,2,M,L,param);
% animateWoidTrajectories(xyarray,['woidlino_test_movies/test_periodic_haptotaxis_' num2str(param.f_hapt)],L);

% angle noise for multiple rods
% test angle noise 
param.angleNoise = 0.02; % not much point making this any bigger than 10, because it's angular
param.k_theta = 0;
xyarray = runWoids(20,40,M,L,param);
filename = ['woidlino_test_movies/test_40rods_' ...
    'angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)];

% % test haptotaxis for multiple rods
% rng(2)
% param.f_hapt = 0.1;
% param.k_theta = 2;
% xyarray = runWoids(50,40,M,L,param);
% animateWoidTrajectories(xyarray,['woidlino_test_movies/test_40rods_haptotaxis_' num2str(param.f_hapt) '_ktheta_' num2str(param.k_theta)],L);

% % test haptotaxis for multiple rods with angle noise
% rng(2)
% param.f_hapt = 0.25;
% param.k_theta = 2;
% param.angleNoise = 1;
% xyarray = runWoids(50,40,M,L,param);
% animateWoidTrajectories(xyarray,['woidlino_test_movies/test_40rods_haptotaxis_' ...
%     num2str(param.f_hapt) '_angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)],L);

% L = [15 15];
% N = 50;
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_noflux_square',L,rc);
% 
% param.bc = 'periodic';
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_periodic_square',L,rc);
% 
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
%     num2str(param.eps_LJ,'%1.0e') ...
%     '_LJhead'],L,param.rc);

% % test long woidlinos without volume exclusion

% % test roaming state
% param.slowingMode = 'gradual';
% param.k_roam = 0.01;
% param.k_unroam = 0.05;
% param.vs = 0.01;
% param.revRateClusterEdge = 1.6;
% xyarray = runWoids(250,N,M,L,param);
% animateWoidTrajectories(xyarray,...
%     ['woidlino_test_movies/test_longBody_periodic_square'...
%     '_eps_LJ_' num2str(param.eps_LJ,'%1.0e'),...
%     '_noVolExcl' '_slowing' param.slowingMode ...
%     '_roam_' num2str(param.k_roam) '_' num2str(param.k_unroam)],L,rc0);

% % test stochastic slowing
% param.dT = rc/param.v0/16; % dT: time step, gets adapted in simulation
% param.saveEvery = round(1/2/param.dT);
% param.slowingMode = 'stochastic_bynode';
% param.k_dwell = 0.0036;%1/275;%1/4;
% param.k_undwell = 1.1;%1/0.9; %1/2.2;
% param.dkdN_dwell = 2;
% param.dkdN_undwell = param.dkdN_dwell;
% param.revRateClusterEdge = 1.6;
% param.vs = 0.014;
% xyarray = runWoids(150,N,M,L,param);
% animateWoidTrajectories(xyarray,...
%     ['woidlino_test_movies/test_longBody_periodic_square'...
%     '_noVolExcl' '_slowing' param.slowingMode ...
%     '_revRateClusterEdge_' num2str(param.revRateClusterEdge) ...
%     '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
%     '_dkdN_' num2str(param.dkdN_dwell)],L,rc0);

% test feeding
L = [7.5 7.5];
T = 7200;
M = 18;
param.rc = 0;
param.r_feed = 1/100;
param.k_unroam = 10;
param.slowingMode = 'stochastic_bynode';
param.k_dwell = 0.0036;
param.k_undwell = 1.1;
param.reversalMode = 'density';
param.revRateClusterEdge = 0;
param.drdN_rev = 0.4;
param.vs = 0.018;
param.dkdN_dwell = 0;
param.dkdN_undwell = 1.4;
param.angleNoise = 0.02;
[xyarray, ~, food] = runWoids(T,40,M,L,'bc','periodic',param);
filename = ['woidlino_test_movies/40rodsM' num2str(M) ...
    '_sweeping_feedrate_' num2str(param.r_feed) '_kunroam_' num2str(param.k_unroam)...
    '_angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)];

% % N2-like
% param.dT = rc/param.v0/16; % dT: time step, gets adapted in simulation
% param.saveEvery = round(1/2/param.dT);
% param.slowingMode = 'stochastic';
% param.k_dwell = 1/4.4;%1/275;%1/4;
% param.k_undwell = 1/2.2;%1/0.9; %1/2.2;
% param.dkdN_dwell = 0;
% param.revRateClusterEdge = 0.25;
% param.v0 = 0.14;
% param.vs = 0.014;
% xyarray = runWoids(150,N,M,L,param);
% animateWoidTrajectories(xyarray,...
%     ['woidlino_test_movies/test_longBody_periodic_square_N2'...
%     '_noVolExcl' '_slowing' param.slowingMode ...
%     '_revRateClusterEdge_' num2str(param.revRateClusterEdge) ...
%     '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
%     '_dkdN_' num2str(param.dkdN_dwell)],L,rc0);

% N = 1;
% param.bc = 'free';
% param.k_dwell = 2;%1/275;%1/4;
% param.k_undwell = 0.1;%1/0.9; %1/2.2;
% param.dkdN_dwell = 0;
% param.revRate = 0.8;
% param.dT = min(1/2,rc0/param.v0/8);
% param.saveEvery = round(1/4/param.dT);
% xyarray = runWoids(50,N,M,L,param);
% animateWoidTrajectories(xyarray,...
%     ['woidlino_test_movies/test_single_periodic_square'...
%     '_noVolExcl' '_slowing' param.slowingMode ...
%     '_revRate_' num2str(param.revRateClusterEdge) ...
%     '_dwell_' num2str(param.k_dwell) '_'
%     num2str(param.k_undwell)],L,rc0);\
%% make movie and other plots
animateWoidTrajectories(xyarray,filename,L,0.035,food);

% pcf_mean = inf_pcf(xyarray,'complexsim',min(param.dT*param.saveEvery/3,1));
% figure
% plot((0.1:0.1:2) - 0.1/2,pcf_mean,'LineWidth',2)
% xlabel('r (mm)'), ylabel('pcf')
% set(gcf,'PaperUnits','centimeters')
% exportfig(gcf,[filename '.eps']);
% system(['epstopdf ' filename '.eps']);
% system(['rm ' filename '.eps']);
% 
% save(['../results/woidlinos/tests/' strrep(filename,'woidlino_test_movies/','') '.mat'])