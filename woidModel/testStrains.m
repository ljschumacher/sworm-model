% test woid model for different strain parameterisations
% to do:

clear

N = 40;
M = 50;
L = 8/2;
dT = 1.2/17/0.33/4/2;
T = 10000;

xyphiarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',pi*37/180);
save('results/DA609_noflux')
animateWoidTrajectories(xyphiarray,'tests/DA609_noflux',L);

xyphiarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',pi*37/180,'slowingNodes',1:M);
save('results/DA609_noflux_slowingNodesAll')
animateWoidTrajectories(xyphiarray,'tests/DA609_noflux_slowingNodesAll',L);

xyphiarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'v0',0.14,...
    'omega_m',2*pi*0.26,...
    'theta_0',pi*37/180,...
    'revRate',1/11,...
    'revTime',1.5,...
    'revRateCluster',1/11,...
    'headNodes',[],...
    'tailNodes',[]);
save('results/N2_noflux')
animateWoidTrajectories(xyphiarray,'tests/N2_noflux',L);




