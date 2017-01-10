% test SPP model

clear
close all

% general model parameters for all test - unless set otherwise
N = 40;
M = 18;
L = 12/2;
dT = 1.2/17/0.33/4;
% single worm
xyphiarray = runWoids(500,1,M,2,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux');

xyphiarray = runWoids(500,1,M,2,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_ktheta0');

xyphiarray = runWoids(500,1,M,2,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_revRate0');

xyphiarray = runWoids(500,1,M,2,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_undulations0');

% two worms
xyphiarray = runWoids(500,2,M,2,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux');

xyphiarray = runWoids(500,2,M,2,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_ktheta0');

xyphiarray = runWoids(500,2,M,2,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_revRate0');

xyphiarray = runWoids(500,2,M,2,'bc','noflux','dT',dT,'rs',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_rslow0');

xyphiarray = runWoids(500,2,M,2,'bc','noflux','dT',dT,'slowingNodes',1:M);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_slowingNodesAll');

xyphiarray = runWoids(500,2,M,2,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux_undulations0');

% many worms
xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux');

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_ktheta0');

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'revRate',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_revRate0');

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'rs',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_rslow0');

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'slowingNodes',1:M);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_slowingNodesAll');

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',0,'omega_m',0,'deltaPhase',0);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux_undulations0');



