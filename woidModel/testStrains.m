% test woid model for different strain parameterisations
% to do:

clear

N = 40;
M = 50;
L = 8/2;
dT = 1.2/M/0.33/4;
T = 10000;
saveevery = 16;

xyphiarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',pi*37/180);
xyphiarray = xyphiarray(:,:,:,1:saveevery:end);
save('results/DA609_noflux')
animateWoidTrajectories(xyphiarray,'tests/DA609_noflux',L);

xyphiarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',pi*37/180,'slowingNodes',1:M);
xyphiarray = xyphiarray(:,:,:,1:saveevery:end);
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
xyphiarray = xyphiarray(:,:,:,1:saveevery:end);
save('results/N2_noflux')
animateWoidTrajectories(xyphiarray,'tests/N2_noflux',L);




