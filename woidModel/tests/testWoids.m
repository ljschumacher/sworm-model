% test SPP model

% issues/todo:
% - fix tests so that straight worm bends and bend worm straightens
% - always reseed random number generator before each simulation?
% - turn this into individual functions? can keep oin same file

clear
close all
addpath('../')
addpath('../visualisation')
addpath('../analysis')
% general model parameters for all test - unless set otherwise
N = 40;
M = 36;
param.dT = 0.035/0.33/16; % baseline timestep, eg rc/v0/8 when bending
% set this to  rc/v0/16 to get better reversal accuracy
param.saveEvery = 32;
food = [];
% single worm
L = [2, 2];
rng(1)

% xyarray = runWoids(10,1,36,L,'bc','periodic',param);
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWormM' num2str(M)],L);
% % plot distribution of lengths to check length conservation
% figure, histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
%     'Normalization','Probability','EdgeColor','none')
% xlabel('L'), ylabel('P')
% set(gcf,'PaperUnits','centimeters')
% movfilename = ['lengthDistribution1Worm'];
% exportfig(gcf,[movfilename '.eps']);
% system(['epstopdf ' movfilename '.eps']);
% system(['rm ' movfilename '.eps']);

% % test N2-like worm
% xyarray = runWoids(15,1,M,L,'bc','free',param,...
%     'v0',0.14,'omega_m',2*pi*0.25);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_N2like',L);

% xyarray = runWoids(20,1,M,L,'bc','noflux',param,'k_theta',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_noflux_ktheta0',L);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux',param,'revRate',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_noflux_revRate0',L);
% 
% rng(6)
% xyarray = runWoids(12,1,M,L,'bc','periodic',param,...
%     'revRate',0.5);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_periodic_reversals',L);

% rng(6)
% tic
% xyarray = runWoids(20,1,M,L,'bc','periodic',param,...
%     'revRate',0.5,'theta_0',0,'omega_m',0,'deltaPhase',0);
% toc
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_periodic_undulations0',L);


% xyarray = runWoids(20,1,M,L,'bc','free',param,...
%     'vs',0.014,'slowingMode','stochastic','k_dwell',1/4,'k_undwell',1/2.2);
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWormM' num2str(M) '_dwelling'],L);
% 
% angleNoise = 1;
% xyarray = runWoids(40,1,M,L,'bc','free',param,...
%     'angleNoise',angleNoise);
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWorm_free_angleNoise' num2str(angleNoise)],L);
% 
% rng(1)
% angleNoise = 1;
% xyarray = runWoids(40,1,18,L,'bc','free',param,...
%     'angleNoise',angleNoise);
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWorm_M18_free_angleNoise' num2str(angleNoise)],L);
%
% xyarray = runWoids(20,1,M,L,'bc','noflux',param,'v0',1e-4,'vs',1e-4,'omega_m',0,'revRate',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_targetcurvatureTest',L,0);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux',param,'v0',1e-4,'vs',1e-4,'theta_0',0,'omega_m',0,'deltaPhase',0,'revRate',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_straighteningTest',L,0);
% 
% % test resumable simulations
% [xyarray1, currentState] = runWoids(5,1,M,L,'bc','noflux',param);
% xyarray2 = runWoids(5,1,M,L,'bc','noflux',param,'resumeState',currentState);
% xyarray = cat(4,xyarray1,xyarray2(:,:,:,2:end));
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWormM' num2str(M) '_noflux_resumed'],L);
% 
% % test feeding
% L = [2.5 2.5];
% [xyarray, ~, food] = runWoids(10,1,M,L,'bc','periodic',param,'r_feed',5);
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWormM' num2str(M) '_feeding'],L,0.035,food);

% % two worms
% L = [2 2];
% xyarray = runWoids(20,2,M,L,'bc','noflux',param);
% animateWoidTrajectories(xyarray,['woid_test_movies/twoWormsM' num2str(M) '_noflux'],L);
% 
% % plot distribution of lengths to check length conservation - doens't work
% % with periodic boundary conditions
% figure, histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
%     'Normalization','Probability','EdgeColor','none')
% xlabel('L'), ylabel('P')
% set(gcf,'PaperUnits','centimeters')
% movfilename = ['lengthDistribution2Worms'];
% exportfig(gcf,[movfilename '.eps']);
% system(['epstopdf ' movfilename '.eps']);
% system(['rm ' movfilename '.eps']);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux',param,'k_theta',00);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_ktheta0',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux',param,'revRate',0,'revRateCluster',0,'revRateClusterEdge',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_revRate0',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux',param,'headNodes',[],'tailNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_revUnresponsive',L);
% 
% rng(1)
% xyarray = runWoids(20,2,M,L,'bc','noflux',param,'revRateClusterEdge',10/13);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_revClusterEdgeIncreased',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux',param,'rc',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_rc0',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux',param,'slowingNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_slowingNodesNone',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux',param,'slowingMode','abrupt');
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_slowingAbrupt',L);
% 
% rng(1)
% xyarray = runWoids(20,3,M,L,'bc','periodic',param,...
%     'vs',1e-2,'slowingMode','abrupt','k_roam',0.1,'k_unroam',0.1);
% animateWoidTrajectories(xyarray,['woid_test_movies/twoWormsM' num2str(M) '_roaming'],L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux',param,...
%     'theta_0',0,'omega_m',0,'deltaPhase',0,'revRate',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_undulations0',L);
% 
% rng(1)
% eps_LJ = 2e-3;
% xyarray = runWoids(40,2,M,L,'bc','periodic',param,...
%     'r_LJcutoff',5*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,'LJnodes',1);
% animateWoidTrajectories(xyarray,['woid_test_movies/twoWorms_periodic_LennardJones' num2str(eps_LJ,'%1.0e') '_head'],L);
%
% rng(2)
% xyarray = runWoids(20,2,M,L,'bc','periodic',param,'f_hapt',1);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_periodic_haptotaxis',L);

% % test density-dependent reversals
% rng(2)
% xyarray = runWoids(20,2,M,L,'bc','periodic',param,...
%     'reversalMode','density','drdN_rev',1);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_periodic_rev_density',L);

% many worms
L = [7.5, 7.5];
rng(1)
% xyarray = runWoids(20,N,M,L,'bc','noflux',param);
% animateWoidTrajectories(xyarray,['woid_test_movies/40wormsM' num2str(M) '_noflux'],L);
% 
% % plot distribution of lengths to check length conservation
% figure, histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
%     'Normalization','Probability','EdgeColor','none')
% xlabel('L'), ylabel('P')
% set(gcf,'PaperUnits','centimeters')
% movfilename = ['lengthDistribution40Worms'];
% exportfig(gcf,[movfilename '.eps']);
% system(['epstopdf ' movfilename '.eps']);
% system(['rm ' movfilename '.eps']);

% xyarray = runWoids(80,N,M,L,'bc','noflux',param,'k_theta',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_ktheta0',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux',param,'revRate',0,'revRateCluster',0,'revRateClusterEdge',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_revRate0',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux',param,'tailNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_revHeadOnly',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux',param,'headNodes',[],'tailNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_revUnresponsive',L);

% % test density-dependent reversals
% xyarray = runWoids(20,40,M,L,'bc','periodic',param,...
%     'reversalMode','density','drdN_rev',0.5);
% animateWoidTrajectories(xyarray,'woid_test_movies/40Worms_periodic_rev_density',L);

% k_l = 54;
% xyarray = runWoids(80,N,M,L,'bc','free','k_l',k_l,param,'revRateClusterEdge',10/13);
% animateWoidTrajectories(xyarray,['woid_test_movies/40Worms_noflux_revClusterEdgeIncreased' '_kl_' num2str(k_l)],L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux',param,'rc',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_rc0',L);

% xyarray = runWoids(80,N,M,L,'bc','noflux',param,'slowingNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_slowingNodesNone',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux',param,...
%     'theta_0',0,'omega_m',0,'deltaPhase',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_undulations0',L);
% 
%

% % % test feeding without volume exclusion
% M = 18;
% param.rc = 0;
% param.k_l = 80;
% param.r_feed = 1/40;
% param.k_unroam = 10;
% param.slowingMode = 'stochastic_bynode';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.reversalMode = 'density';
% param.revRateClusterEdge = 0;
% param.drdN_rev = 0.4;
% param.vs = 0.018;
% param.dkdN_dwell = 0;
% param.dkdN_undwell = 1.4;
% param.angleNoise = 1;
% [xyarray, ~, food] = runWoids(1000,40,M,L,'bc','periodic',param);
% movfilename = ['woid_test_movies/40WormM' num2str(M) ...
%     '_sweeping_feedrate_' num2str(param.r_feed) '_kunroam_' num2str(param.k_unroam)...
%     '_angleNoise' num2str(param.angleNoise)];

% % test feeding without volume exclusion, with haptotaxis
% M = 18;
% param.rc = 0;
% param.k_l = 80;
% param.r_feed = 1/40;
% param.k_unroam = 10;
% param.slowingMode = 'stochastic_bynode';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.revRateClusterEdge = 1;
% param.vs = 0.018;
% param.dkdN_dwell = 0.6;
% param.dkdN_undwell = param.dkdN_dwell;
% param.angleNoise = 1;
% param.f_hapt = 0.2;
% [xyarray, ~, food] = runWoids(1000,40,M,L,'bc','periodic',param);
% animateWoidTrajectories(xyarray,['woid_test_movies/40WormM' num2str(M) ...
%     '_sweeping_feedrate_' num2str(param.r_feed) '_kunroam_' num2str(param.k_unroam)...
%     '_angleNoise' num2str(param.angleNoise) '_haptotaxis_' num2str(param.f_hapt)],L,0.035,food);

% % test feeding without volume exclusion for many worms
% L = 2*[7.5 7.5];
% M = 18;
% param.rc = 0;
% param.k_l = 80;
% param.r_feed = 1/40;
% param.k_unroam = 10;
% param.slowingMode = 'stochastic_bynode';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.revRateClusterEdge = 1;
% param.vs = 0.018;
% param.dkdN_dwell = 0.6;
% param.dkdN_undwell = param.dkdN_dwell;
% param.angleNoise = 1;
% [xyarray, ~, food] = runWoids(1000,200,M,L,'bc','periodic',param);
% animateWoidTrajectories(xyarray,['woid_test_movies/200WormM' num2str(M) ...
%     '_sweeping_feedrate_' num2str(param.r_feed) '_kunroam_' num2str(param.k_unroam)...
%     '_angleNoise' num2str(param.angleNoise)],L,0.035,food);

% test feeding
T = 1000%7200;
bc = 'periodic';
param.k_l = 80;
param.r_feed = 0%1/100;
param.k_unroam = 10;
param.slowingMode = 'stochastic_bynode';
param.k_dwell = 0.0036;
param.k_undwell = 1.1;
param.dkdN_dwell = 0.05;
param.dkdN_undwell = 0.4;
param.reversalMode = 'density';
param.revRateClusterEdge = 0;
param.drdN_rev = 0.05;
param.ri = 5*0.035;
param.vs = 0.018;
param.v0 = 0.33;
% param.omega_m = 2*pi*0.25
param.f_hapt = 0%0.05;
param.eps_LJ = 2e-5;
param.LJmode = 'soft';
param.LJnodes = 1:M;
param.sigma_LJ = 2*0.035;
param.r_LJcutoff = 34*0.035;
[xyarray, currentState, food] = runWoids(T,40,M,L,'bc',bc,param);
movfilename = ['woid_test_movies/40WormM' num2str(M) '_sweeping_feedrate_' num2str(param.r_feed) ...
    '_kunroam_' num2str(param.k_unroam) '_haptotaxis_' num2str(param.f_hapt) ...
    '_LJsoft_' num2str(param.eps_LJ) '_longRange_headOnly'
    '_ri_' num2str(param.ri) ...
    ];

% % test attraction on head-only
% eps_LJ = 5e-3;
% xyarray = runWoids(50,40,M,L,'bc','periodic',param,...
%     'r_LJcutoff',3.75*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,'LJnodes',1,'LJmode','soft',...
%     'slowingNodes',[],'revRateClusterEdge',0,'vs',0.33,'k_l',80 ...
%     );
% animateWoidTrajectories(xyarray,['woid_test_movies/40Worms_periodic_LennardJones' num2str(eps_LJ,'%1.0e')...
%     '_head' '_slowingNodesNone' '_noRev' ...
%     ],L);

% % test volume exclusion only through LJ force
% eps_LJ = 1e-2;
% xyarray = runWoids(20,40,M,L,'bc','periodic',param,...
%     'rc',0,'r_LJcutoff',2*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,'LJnodes',1:M,...
%     'slowingNodes',[],...
%     'revRateClusterEdge',0 ...
%     );
% animateWoidTrajectories(xyarray,['woid_test_movies/40Worms_periodic_LennardJones' num2str(eps_LJ,'%1.0e')...
%     '_asVolExcl'...
%     ],L);

% % test soft lennard-jones potential
% eps_LJ = 1e-1;
% k_l = 320;
% k_theta = 80;
% xyarray = runWoids(500,40,M,L,'bc','periodic',param,...
%     'slowingMode','stochastic_bynode','k_dwell',0,'k_undwell',0,...
%     'revRateClusterEdge',0,'dkdN_dwell',0,...
%     'vs',0.33,'k_l',k_l,'k_theta',k_theta,'r_LJcutoff',3.75*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,...
%     'LJmode','soft'...
%     );
% movfilename = ['woid_test_movies/40Worms_periodic_LennardJones' num2str(eps_LJ,'%1.0e')...
%     '_soft' '_kl' num2str(k_l) '_ktheta' num2str(k_theta)];

% % test soft lennard-jones potential - slow movement
% eps_LJ = 2e-2;
% k_l = 80;
% k_theta = 20;
% v = 0.033;
% xyarray = runWoids(500,40,M,L,'bc','periodic',param,...
%     'slowingMode','stochastic_bynode','k_dwell',0,'k_undwell',0,...
%     'revRateClusterEdge',0,'dkdN_dwell',0,...
%     'v0',v,'omega_m',2*pi*0.6*v/0.33,'vs',v,'k_l',k_l,'k_theta',k_theta,'r_LJcutoff',3.75*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,...
%     'LJmode','soft'...
%     );
% movfilename = ['woid_test_movies/40Worms_periodic_LennardJones' num2str(eps_LJ,'%1.0e')...
%     '_soft' '_kl' num2str(k_l) '_ktheta' num2str(k_theta) '_v' num2str(v)];

% % test soft lennard-jones potential as volume exclusion
% eps_LJ = 2e-2;
% xyarray = runWoids(40,39,M,L,'bc','periodic',param,...
%     'rc',0,'r_LJcutoff',2*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,'LJnodes',1:M,...
%     'LJmode','soft','slowingNodes',[],...
%     'revRateClusterEdge',0 ...
%     );
% movfilename = ['woid_test_movies/40Worms_periodic_LennardJones' num2str(eps_LJ,'%1.0e')...
%     '_soft_asVolExcl'];

% % test stochastic slowing
% param.slowingMode = 'stochastic';
% param.k_dwell = 1/4;
% param.k_undwell = 1/2.2;
% param.vs = 0.014;
% param.bc = 'periodic';
% xyarray = runWoids(20,N,M,L,param);
% movfilename = ['woid_test_movies/40Worms_periodic_square'...
%     '_slowing' param.slowingMode ...
%     '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell)];

% param.bc = 'periodic';
% param.slowingMode = 'stochastic';
% param.k_dwell = 0.0036;%1/275;%1/4;
% param.k_undwell = 1.1;%1/0.9; %1/2.2;
% param.dkdN_dwell = 0.25;
% param.dkdN_undwell = param.dkdN_dwell;
% param.revRateClusterEdge = 1.6;
% param.vs = 0.018;
% xyarray = runWoids(20,N,M,L,param);
% movfilename = ['woid_test_movies/40worms_periodic_square'...
%     '_revRateClusterEdge_' num2str(param.revRateClusterEdge) ...
%     '_slowing' param.slowingMode ...
%     '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
%     '_dkdN_' num2str(param.dkdN_dwell) ];

% xyarray = runWoids(80,N,M,L,'bc','periodic',param);
% movfilename = 'woid_test_movies/40worms_periodic';

% % test haptotaxis
% k_l = 80;
% f_hapt = 0.25;
% xyarray = runWoids(500,40,M,L,'bc','periodic',param,...
%     'f_hapt',f_hapt,'ri',3.75*0.035,'k_l',k_l,'slowingNodes',[],'vs',0.33,'revRateClusterEdge',0);
% movfilename = ['woid_test_movies/40Worms_periodic_haptotaxis_' num2str(f_hapt) '_kl_' num2str(k_l)];

% % test haptotaxis with soft LJ as volume exclusion
% k_l = 80;
% f_hapt = 1;
% eps_LJ = 2e-2;
% xyarray = runWoids(500,40,M,L,'bc','periodic',param,...
%     'r_LJcutoff',2*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,'LJmode','soft',...
%     'f_hapt',f_hapt,'ri',3.75*0.035,'k_l',k_l,'slowingNodes',[],'vs',0.33,'revRateClusterEdge',0);
% movfilename = ['woid_test_movies/40Worms_periodic_haptotaxis_' num2str(f_hapt) '_kl_' num2str(k_l)...
%             '_LJ' num2str(eps_LJ,'%1.0e') '_soft_asVolExcl'];

% % test alignment
% k_l = 80;
% f_align = 1;
% ri = 4*0.035;
% xyarray = runWoids(500,40,M,L,'bc','periodic',param,...
%     'f_align',f_align,'ri',ri,'k_l',k_l,'slowingNodes',[],'vs',0.33,'revRateClusterEdge',0);
% movfilename = ['woid_test_movies/40Worms_periodic_align_' num2str(f_align) ...
%      '_ri_' num2str(ri,2) '_kl_' num2str(k_l)];

% % test aggregation with circular boundary
% L = 4.25;
% param.k_l = 80;
% param.slowingMode = 'stochastic_bynode';
% param.k_dwell = 0.0036;
% param.k_undwell = 1.1;
% param.reversalMode = 'density';
% param.drdN_rev = 0.05;
% param.ri = 5*0.035;
% param.revRateClusterEdge = 0;
% param.vs = 0.018;
% param.dkdN_dwell = 0.05;
% param.dkdN_undwell = 0.4;
% 
% [xyarray, currentState, food] = runWoids(500,40,M,L,'bc','noflux',param);
% movfilename = ['woid_test_movies/40WormM' num2str(M) '_circBound' '_midRange'];

%% make movie and other plots
animateWoidTrajectories(xyarray,movfilename,L,0.035,food);

pcf_mean = inf_pcf(xyarray,'complexsim',min(param.dT*param.saveEvery/3,1));
figure
plot((0.1:0.1:2) - 0.1/2,pcf_mean,'LineWidth',2)
xlabel('r (mm)'), ylabel('pcf')
set(gcf,'PaperUnits','centimeters')
exportfig(gcf,[movfilename '.eps']);
system(['epstopdf ' movfilename '.eps']);
system(['rm ' movfilename '.eps']);

save(['../results/woids/tests/' strrep(movfilename,'woid_test_movies/','') '.mat'])