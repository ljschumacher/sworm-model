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
T = 50; % T: simulation duration (number of time-steps)
param.dT = rc/param.v0/16; % dT: time step, gets adapted in simulation
param.saveEvery = round(1/param.dT);
param.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
% undulations
param.omega_m = 0; % angular frequency of oscillation of movement direction, default 0.6 Hz
param.theta_0 = 0; % amplitude of oscillation of movement direction, default pi/4
param.deltaPhase = 0; % for phase shift in undulations and initial positions, default 0.11
% -- reversal parameters --
param.ri = 3*rc;% ri: radius at which rods register contact (default 3 rc)
% -- slow-down parameters --
param.vs = param.v0;% vs: speed when slowed down (default v0/3)
param.slowingNodes = 1:M;% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
% -- Lennard-Jones parameters --
param.r_LJcutoff = 4*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
param.sigma_LJ = 0;
param.eps_LJ = 0;% eps_LJ: strength of LJ-potential
food = [];

% % test paired initial conditions
% for L = [2.4]
% for rngCtr = 5
% % xyarray2 = runWoids(5,1,M,L,'bc','noflux',param,'resumeState',currentState);
% % xyarray = cat(4,xyarray1,xyarray2(:,:,:,2:end));
% param.bc = 'free';
% param.sigma_LJ = 0;
% param.k_theta = 0;
% param.angleNoise = 0.05;
% param.dT = rc/param.v0/8; % dT: time step, gets adapted in simulation
% param.saveEvery = round(1/param.dT/4);
% param.slowingMode = 'stochastic_bynode';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.dkdN_dwell = 1;
% param.dkdN_undwell = 0.05;
% param.reversalMode = 'density';
% param.drdN_rev = 1;
% param.revRateClusterEdge = 0;
% param.r_LJcutoff = -1
% % set up initial conditions
% rng(rngCtr)
% [~, currentState] = runWoids(1,2,M,[L L],param);
% % continue with random seed
% rng('shuffle')
% xyarray = runWoids(60,2,M,[L L],param,'resumeState',currentState);
% filename = ['woidlino_test_movies/test_paired_' num2str(rngCtr) '_L_' num2str(L)];
% animateWoidTrajectories(xyarray,filename,[L L],0.035,food);
% end
% end
% 
% % test paired initial conditions with haptotaxis
% L = [2.4]
% rngCtr = 5
% % xyarray2 = runWoids(5,1,M,L,'bc','noflux',param,'resumeState',currentState);
% % xyarray = cat(4,xyarray1,xyarray2(:,:,:,2:end));
% param.bc = 'free';
% param.sigma_LJ = 0;
% param.k_theta = 0;
% param.angleNoise = 0.05;
% param.dT = rc/param.v0/8; % dT: time step, gets adapted in simulation
% param.saveEvery = round(1/param.dT/4);
% param.slowingMode = 'stochastic_bynode';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.dkdN_dwell = 0;
% param.dkdN_undwell = 0;
% param.reversalMode = 'density';
% param.drdN_rev = 0;
% param.revRateClusterEdge = 0;
% param.Rif = 1.2/0.035;
% param.f_hapt = 0.2;
% param.haptotaxisMode = 'weighted_additive';
% param.r_LJcutoff = -1
% % set up initial conditions
% rng(rngCtr)
% [~, currentState] = runWoids(1,2,M,[L L],param);
% % continue with random seed
% rng('shuffle')
% xyarray = runWoids(60,2,M,[L L],param,'resumeState',currentState);
% filename = ['woidlino_test_movies/test_paired_' num2str(rngCtr) '_L_' num2str(L)...
%     '_haptotaxis_' param.haptotaxisMode '_' num2str(param.f_hapt) ];
% animateWoidTrajectories(xyarray,filename,[L L],0.035,food);
% 
% % test clustered initial conditions
% rng(1)
% param.bc = 'free';
% param.sigma_LJ = 0;
% L = 1.8;
% param.k_theta = 0;
% param.angleNoise = 0.05;
% param.dT = rc/param.v0/8; % dT: time step, gets adapted in simulation
% param.saveEvery = round(1/param.dT/4);
% param.slowingMode = 'stochastic_bynode';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.dkdN_dwell = 0.05;
% param.dkdN_undwell = 0.05;
% param.reversalMode = 'density';
% param.drdN_rev = 0.05;
% param.revRateClusterEdge = 0;
% param.r_LJcutoff = -1
% xyarray = runWoids(3,40,18,L,param);
% filename = ['woidlino_test_movies/test_clustered_L_' num2str(L(1)) '_Rgyr_' num2str(sqrt(sum(var(xyarray(:,1,:,end)))),2)];

% test angle noise 
N = 40;
M = 18;
T = round(100*0.33/0.14);
maxLag = round(30*0.33/0.14);
figure, hold on
noiseLevels = [0.04, 0.05, 0.06]*sqrt(0.14/0.33);
param.ri = 0; % no interaction
param.v0 = 0.14;
param.vs = 0.14;
param.dT = rc/0.33/8; % dT: time step, gets adapted in simulation
param.saveEvery = round(1/2/param.dT);
for noiseLevel = noiseLevels
    param.angleNoise = noiseLevel;% not much point making this any bigger than 10, because it's angular, and even smaller when worm has no stiffness
    param.bc = 'free';
    param.k_theta = 0;
    L = [3 3];
    rng(1)
    xyarray = runWoids(T,N,M,L,param);
    filename = ['woidlino_test_movies/test_free_M' num2str(M)...
        'angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)...
        '_dT' num2str(param.dT) '_N2'];
    % plot velocity autocorrelation of worm worm trajectory to calibrate
    % persistence
    for nn = 1:N
    velocities = gradient(squeeze(xyarray(nn,1,:,:)))./(param.dT*param.saveEvery); % head orientation
    vac(nn,:) = vectorAutoCorrelation(velocities,round(maxLag./(param.dT*param.saveEvery)));
    end
    vac = mean(vac);
    plot(linspace(0,maxLag,length(vac)),vac)
    % plot average heading change vs time ?
end
ylabel('normalised velocity autocorrelation')
xlabel('time (s)')
h = refline(0,0.23);
h.LineStyle = '--';
h.Color = [0 0 0];
ylim([0 1])
xlim([0 maxLag])
lh = legend(num2str(noiseLevels'),'Location','SouthWest');
lh.Title.String = 'noise \eta';
set(gcf,'PaperUnits','centimeters')
exportfig(gcf,[filename '_vac.eps'],'Color','rgb');
system(['epstopdf ' filename '_vac.eps']);
system(['rm ' filename '_vac.eps']);

% % test haptotaxis
% rng(2)
% L = [3 3];
% param.f_hapt = 0.1;
% param.k_theta = 0;
% param.angleNoise = 0.04;
% xyarray = runWoids(20,2,M,L,param);
% filename = ['woidlino_test_movies/test_periodic_haptotaxis_' num2str(param.f_hapt) ...
%     '_angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)];
% 
% % % angle noise for multiple rods
% % % test angle noise 
% param.angleNoise = 0.02; % not much point making this any bigger than 10, because it's angular
% param.k_theta = 0;
% xyarray = runWoids(20,40,M,L,param);
% filename = ['woidlino_test_movies/test_40rods_' ...
%     'angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)];
% 
% % % % test haptotaxis for multiple rods
% L = [3.75 3.75];
% rng(2)
% param.f_hapt = 0.1;
% param.k_theta = 0;
% param.angleNoise = 0.04;
% xyarray = runWoids(50,10,M,L,param);
% filename = ['woidlino_test_movies/test_10rods_haptotaxis_' num2str(param.f_hapt) ...
%     '_angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)];
% 
% % % % test haptotaxis for multiple rods with angle noise
% rng(1)
% param.ri = 1.2;
% param.f_hapt = 0.2;
% param.haptotaxisMode = 'weighted_additive';
% param.k_theta = 0;
% param.angleNoise = 0.05;
% param.dT = rc/param.v0/8; % dT: time step, gets adapted in simulation
% param.saveEvery = round(1/param.dT);
% param.slowingMode = 'stochastic_weighted';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.dkdN_dwell = 0.05;
% param.dkdN_undwell = 0.05;
% param.reversalMode = 'density_weighted';
% param.drdN_rev = 0.05;
% param.revRateClusterEdge = 0;
% param.r_LJcutoff = -1
% N = 40;
% M = 18
% L = [7.5, 7.5];
% [xyarray, currentState, food] = runWoids(1000,N,M,L,param);
% filename = ['woidlino_test_movies/test_' num2str(N) 'rods_M' num2str(M)...
%     '_angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta) ...
%     '_' param.reversalMode '_ri_' num2str(param.ri) ...
%     '_haptotaxis_' param.haptotaxisMode '_' num2str(param.f_hapt) ];
% 
% % % % test attraction with rc=0 for multiple rods with angle noise
% rng(1)
% param.ri = 1.2;
% param.r_feed = 1/100;
% param.k_unroam = 10;
% param.f_hapt = 0.2;
% param.haptotaxisMode = 'weighted';
% param.k_theta = 0;
% param.angleNoise = 0.05;
% param.dT = rc/param.v0/8; % dT: time step, gets adapted in simulation
% param.saveEvery = round(1/param.dT);
% param.slowingMode = 'stochastic_weighted';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.dkdN_dwell = 0.05;
% param.dkdN_undwell = 0.05;
% param.reversalMode = 'density_weighted';
% param.drdN_rev = 0.05;
% param.revRateClusterEdge = 0;
% param.eps_LJ = 0;%1.5e-5;
% param.LJmode = 'soft';
% param.LJnodes = 1:M;
% param.sigma_LJ = 2*0.035;
% param.r_LJcutoff = -1%1.2
% N = 40;
% L = [7.5, 7.5];
% M = 18;
% [xyarray, currentState, food] = runWoids(7200,N,M,L,param);
% filename = ['woidlino_test_movies/test_' num2str(N) 'rods_M' num2str(M) ...
%     '_sweeping_feedrate_' num2str(param.r_feed) ...
%     '_angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)...
%     '_' param.reversalMode '_ri_' num2str(param.ri) ...
%     '_haptotaxis_' param.haptotaxisMode '_' num2str(param.f_hapt) ...
%     '_LJ' param.LJmode '_' num2str(param.eps_LJ) ...'_longRange' '_headOnly'...
%     ...'_sigma_' num2str(param.sigma_LJ) ...
%     ];
% 
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

% N = 1;
% param.revRate = 0.1;
% param.revTime = 10;
% param.k_theta = 0;
% param.angleNoise = 0.04*sqrt(2);
% param.bc = 'free';
% param.dT = rc/param.v0/8;
% param.saveEvery = round(1/2/param.dT);
% rng(1)
% xyarray = runWoids(100,N,M,L,param);
% filename = ['woidlino_test_movies/test_reversals_dT' num2str(param.dT) ...
%     '_ktheta' num2str(param.k_theta) '_angleNoise' num2str(param.angleNoise)];

% param.revRate = 0;
% param.revRateClusterEdge = 1;
% param.revTime = 5;
% param.headNodes = 1;
% param.tailNodes = 2;
% N = 50;
% xyarray = runWoids(T,N,M,L,param);
% animateWoidTrajectories(xyarray,'woidlino_test_movies/test_periodic_square_reversalsClusterEdge',L,rc);
% 
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
% 
% % test long woidlinos without volume exclusion
% 
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
% 
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

% % test feeding
% L = [7.5 7.5];
% T = 7200;
% M = 18;
% param.rc = 0;
% param.r_feed = 1/100;
% param.k_unroam = 10;
% param.slowingMode = 'stochastic_bynode';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.reversalMode = 'density';
% param.revRateClusterEdge = 0;
% param.drdN_rev = 0.2;
% param.vs = 0.018;
% param.dkdN_dwell = 0.6;
% param.dkdN_undwell = 2;
% param.angleNoise = 0.02;
% param.k_theta = 0;
% [xyarray, ~, food] = runWoids(T,40,M,L,'bc','periodic',param);
% filename = ['woidlino_test_movies/40rodsM' num2str(M) ...
%     '_sweeping_feedrate_' num2str(param.r_feed) '_kunroam_' num2str(param.k_unroam)...
%     '_angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)];

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
%     num2str(param.k_undwell)],L,rc0);
%% make movie and other plots
animateWoidTrajectories(xyarray,filename,L,0.035,food);
% 
% pcf_mean = inf_pcf(xyarray,'simulation',min(param.dT*param.saveEvery/3,1));
% figure
% semilogy((0.1:0.1:2) - 0.1/2,pcf_mean,'LineWidth',2)
% xlabel('r (mm)'), ylabel('pcf'), ylim([0.1 100]), xlim([0 2])
% set(gcf,'PaperUnits','centimeters')
% exportfig(gcf,[filename '.eps']);
% system(['epstopdf ' filename '.eps']);
% system(['rm ' filename '.eps']);
% 
% save(['../results/woidlinos/tests/' strrep(filename,'woidlino_test_movies/','') '.mat'])