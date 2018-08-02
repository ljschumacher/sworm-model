% test SPP model
% to do:
%   add test cases for both square and unequal dimensions
clear

T = 1000;
N = 40;
M = 17;
L = 12;

xyphiarray = run2DSPC(T,N,M,L,'bc','free');
animateChainTrajectories(xyphiarray,'tests/test_bcfree');

xyphiarray = run2DSPC(T,N,M,L,'bc','noflux');
animateChainTrajectories(xyphiarray,'tests/test_noflux');

% xyphiarray = run2DSPC(T,N,M,L,'bc','periodic');
% animateChainTrajectories(xyphiarray,'tests/test_periodic');



