% test SPP model
% to do:
%   add test cases for both square and unequal dimensions
clear

T = 100;
N = 40;
M = 17;
L = 12;

xyphiarray = runWoids(T,N,M,L,'bc','free');
animateWoidTrajectories(xyphiarray,'tests/test_bcfree');

xyphiarray = runWoids(T,N,M,L,'bc','noflux');
animateWoidTrajectories(xyphiarray,'tests/test_noflux');

% xyphiarray = runWoids(T,N,M,L,'bc','periodic');
% animateWoidTrajectories(xyphiarray,'tests/test_periodic');



