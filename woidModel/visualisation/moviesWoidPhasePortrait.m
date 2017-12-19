% make movies from simulation results...

clear
close all

% general model parameters for all test - unless set otherwise
L = 7.5*[1, 1]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
rc = 0.035;
Nval = 40;
revRatesClusterEdge = [6.4]%fliplr([0, 0.1, 0.2, 0.4, 0.8, 1.6]);
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125]);
slowspeeds = [0.018];
attractionStrengths = [1e-2];
slowingMode = 'stochastic';
LJmode = 'soft';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = 0.125%fliplr([0 1./[8 4 2 1]]);
% num_nbr_max_per_node = 3;
paramCombis = combvec(revRatesClusterEdge,speeds,slowspeeds,attractionStrengths,dkdN_dwell_values);
nParamCombis = size(paramCombis,2);
for paramCtr = 1:nParamCombis % can be parfor but might impair movie quality
    revRateClusterEdge = paramCombis(1,paramCtr);
    speed = paramCombis(2,paramCtr);
    slowspeed = paramCombis(3,paramCtr);
    attractionStrength = paramCombis(4,paramCtr);
    dkdN_dwell = paramCombis(5,paramCtr);
    filename = ['woids_N_' num2str(Nval) '_L_' num2str(L(1)) ...
                ...'_noVolExcl'... '_angleNoise' ...
                '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
                '_dkdN_' num2str(dkdN_dwell)...num2str(num_nbr_max_per_node) ...
                '_epsLJ_' num2str(attractionStrength,'%1.0e') '_' LJmode ...
                '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e')...
                '_run1'];
    if exist(['../results/woids/' filename '.mat'],'file')...
            &&~exist(['../movies/woids/' filename '.mp4'],'file')
        out = load(['../results/woids/' filename '.mat']);
        animateWoidTrajectories(out.xyarray,['../movies/woids/' filename],L,rc);
    elseif ~exist(['../results/woids/' filename '.mat'],'file')
        disp(['no results for ' filename])
    end
end
