% test SPP model

% issues/todo:
% - fix tests so that straight worm bends and bend worm straightens
% - always reseed random number generator before each simulation?

clear
close all

% general model parameters for all test - unless set otherwise
N = 40;
M = 49;
dT = 0.025; % baseline timestep, eg rc/v0/4
saveevery = 4;

% single worm
L = 2;
xyarray = runWoids(50,1,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/singleWorm_noflux',L);
% plot distribution of lengths to check length conservation
histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
    'Normalization','Probability','EdgeColor','none')
xlabel('L'), ylabel('P')
set(gcf,'PaperUnits','centimeters')
filename = ['tests/lengthDistribution1Worm'];
exportfig(gcf,[filename '.eps']);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/singleWorm_noflux_ktheta0',L);

xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/singleWorm_noflux_revRate0',L);

xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/singleWorm_noflux_undulations0',L);

xyarray = runWoids(20,1,M,L,'bc','noflux','dT',dT,'v0',1e-4,'vs',1e-4,'omega_m',0,'revRate',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/singleWorm_targetcurvatureTest',L,0.01);

% two worms
xyarray = runWoids(50,2,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/twoWorms_noflux',L);

% plot distribution of lengths to check length conservation
histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
    'Normalization','Probability','EdgeColor','none')
xlabel('L'), ylabel('P')
set(gcf,'PaperUnits','centimeters')
filename = ['tests/lengthDistribution2Worms'];
exportfig(gcf,[filename '.eps']);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/twoWorms_noflux_ktheta0',L);

xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/twoWorms_noflux_revRate0',L);

xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'headNodes',[],'tailNodes',[]);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/twoWorms_noflux_revUnresponsive',L);

xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'revRateCluster',1/13);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/twoWorms_noflux_revClusterUnreduced',L);

xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'rc',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/twoWorms_noflux_rc0',L);

xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'slowingNodes',[]);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/twoWorms_noflux_slowingNodesNone',L);

xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,'slowingNodes',1:M);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/twoWorms_noflux_slowingNodesAll',L);

xyarray = runWoids(20,2,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0,'revRate',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/twoWorms_noflux_undulations0',L);

% many worms
L = 8.5/2;

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux',L);

% plot distribution of lengths to check length conservation
histogram(squeeze(sum(sqrt(sum(diff(xyarray(:,:,1:2,:),1,2).^2,3)),2)),...
    'Normalization','Probability','EdgeColor','none')
xlabel('L'), ylabel('P')
set(gcf,'PaperUnits','centimeters')
filename = ['tests/lengthDistribution40Worms'];
exportfig(gcf,[filename '.eps']);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux_ktheta0',L);

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux_revRate0',L);

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'tailNodes',[]);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux_revHeadOnly',L);

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'headNodes',[],'tailNodes',[]);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux_revUnresponsive',L);

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'revRateCluster',1/13);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux_revClusterUnreduced',L);

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'rc',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux_rc0',L);

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'slowingNodes',[]);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux_slowingNodesNone',L);

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,'slowingNodes',1:M);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux_slowingNodesAll',L);

xyarray = runWoids(80,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyarray(:,:,:,1:saveevery:end),'tests/40worms_noflux_undulations0',L);
