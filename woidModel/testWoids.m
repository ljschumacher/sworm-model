% test SPP model

% - fix tests so that straight worm bends and bend worm straightens
clear
close all

% general model parameters for all test - unless set otherwise
N = 40;
M = 50;
dT = 1.2/M/0.33/4;

% single worm
L = 2;
xyphiarray = runWoids(500,1,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux',L);
% plot distribution of lengths to check length conservation
histogram(squeeze(sum(sqrt(sum(diff(xyphiarray(:,:,1:2,:),1,2).^2,3)),2)),...
    'Normalization','Probability','EdgeColor','none')
xlabel('L'), ylabel('P')
set(gcf,'PaperUnits','centimeters')
filename = ['tests/lengthDistribution1Worm'];
exportfig(gcf,[filename '.eps']);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

xyphiarray = runWoids(1000,1,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_ktheta0',L);

xyphiarray = runWoids(1000,1,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_revRate0',L);

xyphiarray = runWoids(1000,1,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_undulations0',L);

xyphiarray = runWoids(1000,1,M,L,'bc','noflux','dT',dT,'v0',1e-4,'vs',1e-4,'omega_m',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_targetcurvatureTest',L);

% two worms
xyphiarray = runWoids(1000*2,2,M,L,'bc','noflux','dT',dT/2,'kl',40);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux',L);

% plot distribution of lengths to check length conservation
histogram(squeeze(sum(sqrt(sum(diff(xyphiarray(:,:,1:2,:),1,2).^2,3)),2)),...
    'Normalization','Probability','EdgeColor','none')
xlabel('L'), ylabel('P')
set(gcf,'PaperUnits','centimeters')
filename = ['tests/lengthDistribution2Worms'];
exportfig(gcf,[filename '.eps']);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

xyphiarray = runWoids(1000,2,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_ktheta0',L);

xyphiarray = runWoids(1000,2,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_revRate0',L);

xyphiarray = runWoids(1000,2,M,L,'bc','noflux','dT',dT,'headNodes',[],'tailNodes',[]);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_revUnresponsive',L);

xyphiarray = runWoids(1000,2,M,L,'bc','noflux','dT',dT,'revRateCluster',1/13);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_revClusterUnreduced',L);

xyphiarray = runWoids(1000,2,M,L,'bc','noflux','dT',dT,'rs',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_rslow0',L);

xyphiarray = runWoids(1000,2,M,L,'bc','noflux','dT',dT,'slowingNodes',1:M);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_slowingNodesAll',L);

xyphiarray = runWoids(1000,2,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_undulations0',L);

% many worms
L = 8/2;
xyphiarray = runWoids(4000,N,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux',L);

% plot distribution of lengths to check length conservation
histogram(squeeze(sum(sqrt(sum(diff(xyphiarray(:,:,1:2,:),1,2).^2,3)),2)),...
    'Normalization','Probability','EdgeColor','none')
xlabel('L'), ylabel('P')
set(gcf,'PaperUnits','centimeters')
filename = ['tests/lengthDistribution40Worms'];
exportfig(gcf,[filename '.eps']);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

xyphiarray = runWoids(4000,N,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_ktheta0',L);

xyphiarray = runWoids(4000,N,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_revRate0',L);

xyphiarray = runWoids(4000,N,M,L,'bc','noflux','dT',dT,'headNodes',[],'tailNodes',[]);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_revUnresponsive',L);

xyphiarray = runWoids(4000,N,M,L,'bc','noflux','dT',dT,'revRateCluster',1/13);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_revClusterUnreduced',L);

xyphiarray = runWoids(4000,N,M,L,'bc','noflux','dT',dT,'rs',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_rslow0',L);

xyphiarray = runWoids(4000,N,M,L,'bc','noflux','dT',dT,'slowingNodes',1:M);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_slowingNodesAll',L);

xyphiarray = runWoids(4000,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_undulations0',L);



