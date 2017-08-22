% test SPP model

% issues/todo:
% - fix tests so that straight worm bends and bend worm straightens
% - always reseed random number generator before each simulation?

clear
close all
addpath('visualisation')
% general model parameters for all test - unless set otherwise
N = 40;
M = 49;
dT = 0.035/0.33/16; % baseline timestep, eg rc/v0/8 when bending
% set this to  rc/v0/16 to get better reversal accuracy
saveEvery = 16;

% single worm
L = 2;
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery);
% animateWoidTrajectories(xyarray,'tests/singleWorm_noflux',L);
% % plot distribution of lengths to check length conservation
% histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
%     'Normalization','Probability','EdgeColor','none')
% xlabel('L'), ylabel('P')
% set(gcf,'PaperUnits','centimeters')
% filename = ['tests/lengthDistribution1Worm'];
% exportfig(gcf,[filename '.eps']);
% system(['epstopdf ' filename '.eps']);
% system(['rm ' filename '.eps']);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'k_theta',0);
% animateWoidTrajectories(xyarray,'tests/singleWorm_noflux_ktheta0',L);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRate',0);
% animateWoidTrajectories(xyarray,'tests/singleWorm_noflux_revRate0',L);
% 
% rng(6)
% xyarray = runWoids(12,1,M,[L L],'bc','periodic','dT',dT,'saveEvery',saveEvery,...
%     'revRate',0.5,'theta_0',0,'omega_m',0,'deltaPhase',0);
% animateWoidTrajectories(xyarray,'tests/singleWorm_periodic_undulations0',[L L]);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'v0',1e-4,'vs',1e-4,'omega_m',0,'revRate',0);
% animateWoidTrajectories(xyarray,'tests/singleWorm_targetcurvatureTest',L,0);
% 
% xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'v0',1e-4,'vs',1e-4,'theta_0',0,'omega_m',0,'deltaPhase',0,'revRate',0);
% animateWoidTrajectories(xyarray,'tests/singleWorm_straighteningTest',L,0);
% 
% % two worms
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux',L);
% 
% % plot distribution of lengths to check length conservation
% histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
%     'Normalization','Probability','EdgeColor','none')
% xlabel('L'), ylabel('P')
% set(gcf,'PaperUnits','centimeters')
% filename = ['tests/lengthDistribution2Worms'];
% exportfig(gcf,[filename '.eps']);
% system(['epstopdf ' filename '.eps']);
% system(['rm ' filename '.eps']);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'k_theta',0);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux_ktheta0',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRate',0,'revRateCluster',0,'revRateClusterEdge',0);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux_revRate0',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'headNodes',[],'tailNodes',[]);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux_revUnresponsive',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRateClusterEdge',10/13);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux_revClusterEdgeIncreased',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'rc',0);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux_rc0',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'slowingNodes',[]);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux_slowingNodesNone',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'slowingNodes',1:M);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux_slowingNodesAll',L);
% 
% xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,...
%     'theta_0',0,'omega_m',0,'deltaPhase',0,'revRate',0);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux_undulations0',L);
% 
rng(1)
eps_LJ = 4e-4;
xyarray = runWoids(20,2,M,[L L],'bc','periodic','dT',dT,'saveEvery',saveEvery,...
    'r_LJcutoff',5*0.035,'eps_LJ',eps_LJ,'sigma_LJ',2*0.035);
animateWoidTrajectories(xyarray,['tests/twoWorms_periodic_LennardJones' num2str(eps_LJ,'%1.0e')],[L L]);

% many worms
% L = 8.5/2;
% tic
% rng(1)
% xyarray = runWoids(5,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery);
% toc
% animateWoidTrajectories(xyarray,'tests/40worms_noflux',L);

% % plot distribution of lengths to check length conservation
% histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
%     'Normalization','Probability','EdgeColor','none')
% xlabel('L'), ylabel('P')
% set(gcf,'PaperUnits','centimeters')
% filename = ['tests/lengthDistribution40Worms'];
% exportfig(gcf,[filename '.eps']);
% system(['epstopdf ' filename '.eps']);
% system(['rm ' filename '.eps']);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'k_theta',0);
% animateWoidTrajectories(xyarray,'tests/40worms_noflux_ktheta0',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRate',0,'revRateCluster',0,'revRateClusterEdge',0);
% animateWoidTrajectories(xyarray,'tests/40worms_noflux_revRate0',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'tailNodes',[]);
% animateWoidTrajectories(xyarray,'tests/40worms_noflux_revHeadOnly',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'headNodes',[],'tailNodes',[]);
% animateWoidTrajectories(xyarray,'tests/40worms_noflux_revUnresponsive',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'revRateClusterEdge',10/13);
% animateWoidTrajectories(xyarray,'tests/twoWorms_noflux_revClusterEdgeIncreased',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'rc',0);
% animateWoidTrajectories(xyarray,'tests/40worms_noflux_rc0',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'slowingNodes',[]);
% animateWoidTrajectories(xyarray,'tests/40worms_noflux_slowingNodesNone',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,'slowingNodes',1:M);
% animateWoidTrajectories(xyarray,'tests/40worms_noflux_slowingNodesAll',L);
% 
% xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'saveEvery',saveEvery,...
%     'theta_0',0,'omega_m',0,'deltaPhase',0);
% animateWoidTrajectories(xyarray,'tests/40worms_noflux_undulations0',L);

% L = [7.5, 7.5];
% xyarray = runWoids(80,N,M,L,'bc','periodic','dT',dT,'saveEvery',saveEvery);
% animateWoidTrajectories(xyarray,'tests/40worms_periodic',L);