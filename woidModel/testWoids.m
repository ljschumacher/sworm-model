% test SPP model
% to do:
%   add test cases for both square and unequal dimensions
clear

N = 2;
M = 18;
L = 12;
dT = 1.2/17/0.33/4;

xyphiarray = runWoids(500,N,M,L,'bc','noflux','dT',dT);
animateWoidTrajectories(xyphiarray,'tests/test_noflux');

% xyphiarray = runWoids(500,N,M,L,'bc','free','dT',dT);
% animateWoidTrajectories(xyphiarray,'tests/test_bcfree');

% xyphiarray = runWoids(100,N,M,L,'bc','periodic');
% animateWoidTrajectories(xyphiarray,'tests/test_periodic');



