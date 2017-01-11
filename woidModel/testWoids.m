% test SPP model

clear
close all

% general model parameters for all test - unless set otherwise
N = 40;
M = 18;
dT = 1.2/17/0.33/4;
% single worm
L = 2;
xyphiarray = runWoids(500,1,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux',L);
% plot distribution of lengths to check length conservation
histogram(squeeze(sum(sqrt(sum(diff(xyphiarray(:,:,1:2,:),1,2).^2,3)),2)),...
    'Normalization','Probability','EdgeColor','none')
xlabel('L'), ylabel('P')
set(gcf,'PaperUnits','centimeters')
filename = ['tests/lengthDistribution'];
exportfig(gcf,[filename '.eps']);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

xyphiarray = runWoids(500,1,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_ktheta0',L);

xyphiarray = runWoids(500,1,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_revRate0',L);

xyphiarray = runWoids(500,1,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_undulations0',L);

xyphiarray = runWoids(500,1,M,L,'bc','noflux','dT',dT,'v0',1e-4,'vs',1e-4,'omega_m',1e-4);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_straighteningTest',L);

xyphiarray = runWoids(500,1,M,L,'bc','noflux','dT',dT,'v0',1e-4,'vs',1e-4);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_targetcurvatureTest',L);

% two worms
xyphiarray = runWoids(500,2,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux',L);

xyphiarray = runWoids(500,2,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_ktheta0',L);

xyphiarray = runWoids(500,2,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_revRate0',L);

xyphiarray = runWoids(500,2,M,L,'bc','noflux','dT',dT,'rs',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_rslow0',L);

xyphiarray = runWoids(500,2,M,L,'bc','noflux','dT',dT,'slowingNodes',1:M);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_slowingNodesAll',L);

xyphiarray = runWoids(500,2,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_undulations0',L);

% many worms
L = 12/2;
xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux',L);

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_ktheta0',L);

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_revRate0',L);

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'rs',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_rslow0',L);

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'slowingNodes',1:M);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_slowingNodesAll',L);

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_undulations0',L);



