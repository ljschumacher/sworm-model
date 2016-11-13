% test SPP model
% to do:
%   add test cases for both square and unequal dimensions
clear

N = 40;
M = 17;
L = 12;

xyphiarray = runWoids(500,N,M,L,'bc','free');
animateWoidTrajectories(xyphiarray,'tests/test_bcfree');

xyphiarray = runWoids(1000,N,M,L,'bc','noflux');
animateWoidTrajectories(xyphiarray,'tests/test_noflux');

% xyphiarray = runWoids(100,N,M,L,'bc','periodic');
% animateWoidTrajectories(xyphiarray,'tests/test_periodic');



