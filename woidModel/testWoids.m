% test SPP model
% to do:
%   add test cases for both square and unequal dimensions
clear

N = 40;
M = 16;
L = 12;
dT = 1/20;

xyphiarray = runWoids(250,N,M,L,'bc','free','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/test_bcfree');

xyphiarray = runWoids(1000,N,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/test_noflux');

% xyphiarray = runWoids(100,N,M,L,'bc','periodic');
% animateWoidTrajectories(xyphiarray,'tests/test_periodic');



