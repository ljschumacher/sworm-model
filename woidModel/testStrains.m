% test woid model for different strain parameterisations
% to do:
%   parameterise slow-down
clear

N = 100;
M = 18;
L = 12;
dT = 1.2/17/0.33/4;
T = 3000;

xyphiarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'theta_0',pi*40/180);
animateWoidTrajectories(xyphiarray,'tests/DA609_noflux');

xyphiarray = runWoids(T,N,M,L,'bc','noflux','dT',dT,...
    'v0',0.14,...
    'omega_m',2*pi*0.26,...
    'theta_0',pi*37/180,...
    'deltaPhase',2*pi/M*1.2/0.68,...
    'revRate',1/11,...
    'revTime',1.5);
animateWoidTrajectories(xyphiarray,'tests/N2_noflux');




