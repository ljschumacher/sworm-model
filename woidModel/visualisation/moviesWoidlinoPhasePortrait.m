% make movies from simulation results...

clear
close all

% general model parameters for all test - unless set otherwise
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
% rc = 0.035;
N = 40;
M = 18;
revRatesClusterEdge = 1.6%fliplr([0, 0.4, 0.8, 1.6, 3.2, 6.4]);
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125]);
slowspeeds = [0.018];
attractionStrengths = [0];
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = 1/8%[0 1./[8 4 2 1]];
paramCombis = combvec(revRatesClusterEdge,speeds,slowspeeds,...
    attractionStrengths,dkdN_dwell_values);
nParamCombis = size(paramCombis,2);
angleNoise = 0.02;

for paramCtr = 1:nParamCombis % can be parfor but might impair movie quality
    revRateClusterEdge = paramCombis(1,paramCtr);
    speed = paramCombis(2,paramCtr);
    slowspeed = paramCombis(3,paramCtr);
    attractionStrength = paramCombis(4,paramCtr);
    dkdN_dwell = paramCombis(5,paramCtr);
    filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
        '_noVolExcl' '_angleNoise_' num2str(angleNoise) ...
        '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
        '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
        '_dkdN_' num2str(dkdN_dwell) ...
        '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
        '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e') ...
        '_run1'];
    if exist(['../results/woidlinos/' filename '.mat'],'file')...
            &&~exist(['../movies/woidlinos/' filename '.mp4'],'file')
        out = load(['../results/woidlinos/' filename '.mat']);
        animateWoidTrajectories(out.xyarray,['../movies/woidlinos/' filename],L);%,out.param.rc);
    elseif ~exist(['../results/woidlinos/' filename '.mat'],'file')
        disp(['no results for ' filename])
    end
end
