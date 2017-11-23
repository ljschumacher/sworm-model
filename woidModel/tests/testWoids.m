% test SPP model

% issues/todo:
% - fix tests so that straight worm bends and bend worm straightens
% - always reseed random number generator before each simulation?

clear
close all
addpath('../')
addpath('../visualisation')
% general model parameters for all test - unless set otherwise
N = 40;
M = 49;
dT = 0.035/0.33/16; % baseline timestep, eg rc/v0/8 when bending
% set this to  rc/v0/16 to get better reversal accuracy
saveEvery = 16;

% single worm
L = 2;
rng(1)
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery);
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWormM' num2str(M) '_noflux'],L);
% % plot distribution of lengths to check length conservation
% figure, histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
%     'Normalization','Probability','EdgeColor','none')
% xlabel('L'), ylabel('P')
% set(gcf,'PaperUnits','centimeters')
% filename = ['lengthDistribution1Worm'];
% exportfig(gcf,[filename '.eps']);
% system(['epstopdf ' filename '.eps']);
% system(['rm ' filename '.eps']);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'k_theta',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_noflux_ktheta0',L);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRate',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_noflux_revRate0',L);
% 
% rng(6)
% xyarray = runWoids(12,1,M,[L L],'bc','periodic','dT',dT,'saveEvery',saveEvery,...
%     'revRate',0.5);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_periodic_reversals',[L L]);
%
% rng(6)
% xyarray = runWoids(12,1,M,[L L],'bc','periodic','dT',dT,'saveEvery',saveEvery,...
%     'revRate',0.5,'theta_0',0,'omega_m',0,'deltaPhase',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_periodic_undulations0',[L L]);
%
% xyarray = runWoids(20,1,M,[L L],'bc','free','dT',dT,'saveEvery',saveEvery,...
%     'vs',0.014,'slowingMode','stochastic','k_dwell',1/4,'k_undwell',1/2.2);
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWormM' num2str(M) '_dwelling'],[L L]);
%
% angleNoise = 0.1;
% xyarray = runWoids(40,1,M,[L L],'bc','free','dT',dT,'saveEvery',saveEvery,...
%     'angleNoise',angleNoise);
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWorm_free_angleNoise' num2str(angleNoise)],[L L]);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'v0',1e-4,'vs',1e-4,'omega_m',0,'revRate',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_targetcurvatureTest',L,0);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'v0',1e-4,'vs',1e-4,'theta_0',0,'omega_m',0,'deltaPhase',0,'revRate',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/singleWorm_straighteningTest',L,0);
%
% test resumable simulations
% [xyarray1, currentState] = runWoids(5,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery);
% xyarray2 = runWoids(5,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'resumeState',currentState);
% xyarray = cat(4,xyarray1,xyarray2(:,:,:,2:end));
% animateWoidTrajectories(xyarray,['woid_test_movies/singleWormM' num2str(M) '_noflux_resumed'],L);
%
% two worms
L = [2 2];
% xyarray = runWoids(20,2,M,L,'bc','periodic','dT',dT,'saveEvery',saveEvery);
% animateWoidTrajectories(xyarray,['woid_test_movies/twoWormsM' num2str(M) '_noflux'],L);
% 
% % plot distribution of lengths to check length conservation
% figure, histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
%     'Normalization','Probability','EdgeColor','none')
% xlabel('L'), ylabel('P')
% set(gcf,'PaperUnits','centimeters')
% filename = ['lengthDistribution2Worms'];
% exportfig(gcf,[filename '.eps']);
% system(['epstopdf ' filename '.eps']);
% system(['rm ' filename '.eps']);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'k_theta',00);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_ktheta0',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRate',0,'revRateCluster',0,'revRateClusterEdge',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_revRate0',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'headNodes',[],'tailNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_revUnresponsive',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRateClusterEdge',10/13);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_revClusterEdgeIncreased',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'rc',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_rc0',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'slowingNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_slowingNodesNone',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'slowingMode','abrupt');
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_slowingAbrupt',L);
% 
% rng(1)
% xyarray = runWoids(20,3,M,L,'bc','periodic','dT',dT,'saveEvery',saveEvery,...
%     'vs',1e-2,'slowingMode','abrupt','k_roam',0.1,'k_unroam',0.1);
% animateWoidTrajectories(xyarray,['woid_test_movies/twoWormsM' num2str(M) '_roaming'],L);
%
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,...
%     'theta_0',0,'omega_m',0,'deltaPhase',0,'revRate',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_undulations0',L);
% 
% rng(1)
% eps_LJ = 2e-3;
% xyarray = runWoids(40,2,M,L,'bc','periodic','dT',dT,'saveEvery',saveEvery,...
%     'r_LJcutoff',5*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,'LJnodes',1);
% animateWoidTrajectories(xyarray,['woid_test_movies/twoWorms_periodic_LennardJones' num2str(eps_LJ,'%1.0e') '_head'],L);

% many worms
L = 8.5/2;
rng(1)
% xyarray = runWoids(20,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery);
% animateWoidTrajectories(xyarray,['woid_test_movies/40wormsM' num2str(M) '_noflux'],L);
% 
% % plot distribution of lengths to check length conservation
% figure, histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
%     'Normalization','Probability','EdgeColor','none')
% xlabel('L'), ylabel('P')
% set(gcf,'PaperUnits','centimeters')
% filename = ['lengthDistribution40Worms'];
% exportfig(gcf,[filename '.eps']);
% system(['epstopdf ' filename '.eps']);
% system(['rm ' filename '.eps']);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'k_theta',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_ktheta0',L);

% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRate',0,'revRateCluster',0,'revRateClusterEdge',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_revRate0',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'tailNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_revHeadOnly',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'headNodes',[],'tailNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_revUnresponsive',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRateClusterEdge',10/13);
% animateWoidTrajectories(xyarray,'woid_test_movies/twoWorms_noflux_revClusterEdgeIncreased',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'rc',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_rc0',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'slowingNodes',[]);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_slowingNodesNone',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,...
%     'theta_0',0,'omega_m',0,'deltaPhase',0);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_noflux_undulations0',L);
% 
rng(1)
eps_LJ = 4e-3;
L = [7.5, 7.5];
xyarray = runWoids(100,40,M,L,'bc','periodic','dT',dT,'saveEvery',saveEvery,...
    'r_LJcutoff',5*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,'LJnodes',1,...
    'slowingNodes',[],...
    'revRate', 0, 'revRateCluster', 0,'revRateClusterEdge',0 ...
    ,'theta_0',0,'omega_m',0,'deltaPhase',0 ...
    );
animateWoidTrajectories(xyarray,['woid_test_movies/40Worms_periodic_LennardJones' num2str(eps_LJ,'%1.0e')...
    '_head' '_slowingNodesNone' '_noRev' ...
    ,'_undulations0'...
    ],L);
% 
% rng(1)
% eps_LJ = 1e-2;
% L = [7.5, 7.5];
% xyarray = runWoids(10,40,M,L,'bc','periodic','dT',dT,'saveEvery',saveEvery,...
%     'rc',0,'r_LJcutoff',2*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035,'LJnodes',1:M,...
%     'slowingNodes',[],...
%     'revRateClusterEdge',0 ...
%     );
% animateWoidTrajectories(xyarray,['woid_test_movies/40Worms_periodic_LennardJones' num2str(eps_LJ,'%1.0e')...
%     '_asVolExcl'...
%     ],L);

% test stochastic slowing
% % param.slowingMode = 'stochastic';
% param.k_dwell = 1/4;
% param.k_undwell = 1/2.2;
% param.vs = 0.014;
% xyarray = runWoids(20,N,M,L,param,'saveEvery',saveEvery);
% animateWoidTrajectories(xyarray,...
%     ['woid_test_movies/40Worms_periodic_square'...
%     '_slowing' param.slowingMode ...
%     '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell)],[L L]);

% rng(1)
% L = [7.5, 7.5];
% param.bc = 'periodic';
% param.slowingMode = 'stochastic';
% param.k_dwell = 0.0036;%1/275;%1/4;
% param.k_undwell = 1.1;%1/0.9; %1/2.2;
% param.dkdN_dwell = 0.25;
% param.revRateClusterEdge = 1.6;
% param.vs = 0.018;
% xyarray = runWoids(20,N,M,L,param);
% animateWoidTrajectories(xyarray,...
%     ['woid_test_movies/40worms_periodic_square'...
%     '_revRateClusterEdge_' num2str(param.revRateClusterEdge) ...
%     '_slowing' param.slowingMode ...
%     '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
%     '_dkdN_' num2str(param.dkdN_dwell) ],L);

% L = [7.5, 7.5];
% xyarray = runWoids(80,N,M,L,'bc','periodic','dT',dT,'saveEvery',saveEvery);
% animateWoidTrajectories(xyarray,'woid_test_movies/40worms_periodic',L);