% test woid model for different strain parameterisations
% to do:

clear

N = 40;
M = 49;
L = 8/2;
dT = 1.2/M/0.33/8/2;
T = 10000;
saveevery = 24;

xyarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',pi*37/180);
xyarray = xyarray(:,:,:,1:saveevery:end);
save('results/DA609_noflux')
animateWoidTrajectories(xyarray,'tests/DA609_noflux',L);

xyarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',pi*37/180,'slowingNodes',1:M);
xyarray = xyarray(:,:,:,1:saveevery:end);
save('results/DA609_noflux_slowingNodesAll')
animateWoidTrajectories(xyarray,'tests/DA609_noflux_slowingNodesAll',L);

xyarray = runWoids(4*T,N,M,L,'bc','noflux','dT',dT/4,...
    'theta_0',pi*37/180,'r_LJcutoff',4*0.035,'eps_LJ',2e-6);
xyarray = xyarray(:,:,:,1:saveevery*4:end);
save('results/DA609_noflux_lennardjones2e-6')
animateWoidTrajectories(xyarray,'tests/DA609_noflux_lennardjones2e-6',L);

xyarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'v0',0.14,...
    'omega_m',2*pi*0.26,...
    'theta_0',pi*37/180,...
    'revRate',1/11,...
    'revTime',1.5,...
    'revRateCluster',1/11,...
    'headNodes',[],...
    'tailNodes',[]);
xyarray = xyarray(:,:,:,1:saveevery:end);
save('results/N2_noflux')
animateWoidTrajectories(xyarray,'tests/N2_noflux',L);
