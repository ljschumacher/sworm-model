% test SPP model
% to do:
% - write suite of tests for various cases including different number of
% worms and which interactions are included
clear
close all

% general model parameters for all test - unless set otherwise
N = 40;
M = 18;
L = 12/2;
dT = 1.2/17/0.33/4;

xyphiarray = runWoids(500,1,M,3,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux');

xyphiarray = runWoids(500,1,M,3,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/singleWorm_noflux_ktheta0');

xyphiarray = runWoids(500,2,M,3,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/twoWorms_noflux');

xyphiarray = runWoids(500,2,M,3,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/twoWorm_noflux_ktheta0');

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/40worms_noflux');

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'k_theta',0);
animateWoidTrajectories(xyphiarray,'tests/test_noflux_ktheta0');

xyphiarray = runWoids(3000,N,M,L,'bc','noflux','dT',dT,'rc',0);
animateWoidTrajectories(xyphiarray,'tests/test_noflux_nocontact');



