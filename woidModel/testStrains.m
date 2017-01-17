% test woid model for different strain parameterisations
% to do:
%   - save simulation results
clear

N = 40;
M = 50;
L = 8/2;
dT = 1.2/17/0.33/4/2;
T = 10000;

xyphiarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',pi*37/180);
save('results/DA609_noflux','xyphiarray','L','dT')
animateWoidTrajectories(xyphiarray,'tests/DA609_noflux');

xyphiarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'v0',0.14,...
    'omega_m',2*pi*0.26,...
    'theta_0',pi*37/180,...
    'revRate',1/11,...
    'revTime',1.5,...
    'revRateCluster',1/11,...
    'headNodes',[],...
    'tailNodes',[]);
save('results/N2_noflux','xyphiarray','L','dT')
animateWoidTrajectories(xyphiarray,'tests/N2_noflux');




