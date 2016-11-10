% test SPP model
clear

T = 100;
N = 40;
L = 10;

xyphiarray = run2DSPP(T,N,L,'bc','free');
animateTrajectories(xyphiarray,'tests/test_bcfree');

xyphiarray = run2DSPP(T,N,L,'bc','noflux');
animateTrajectories(xyphiarray,'tests/test_noflux');

xyphiarray = run2DSPP(T,N,L,'bc','periodic');
animateTrajectories(xyphiarray,'tests/test_periodic');



