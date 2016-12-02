% test SPP model
% to do:
%   add test cases for both square and unequal dimensions
clear

N = 40;
M = 18;
L = 12;
dT = 1.2/17/0.33/4;
vs = 0.33/8;

xyphiarray = runWoids(500,N,M,L,'bc','free','dT',dT,'vs',vs);
animateWoidTrajectories(xyphiarray,'tests/test_bcfree');

xyphiarray = runWoids(500,N,M,L,'bc','noflux','dT',dT,'vs',vs);
animateWoidTrajectories(xyphiarray,'tests/test_noflux');

% xyphiarray = runWoids(100,N,M,L,'bc','periodic');
% animateWoidTrajectories(xyphiarray,'tests/test_periodic');



