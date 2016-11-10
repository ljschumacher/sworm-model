% test SPP model
clear

T = 1;
N = 40;
M = 10;
L = 10;

xyphiarray = run2DSPC(T,N,M,L,'bc','free');
animateChainTrajectories(xyphiarray,'tests/test_bcfree');

xyphiarray = run2DSPC(T,N,M,L,'bc','noflux');
animateChainTrajectories(xyphiarray,'tests/test_noflux');

xyphiarray = run2DSPC(T,N,M,L,'bc','periodic');
animateChainTrajectories(xyphiarray,'tests/test_periodic');



