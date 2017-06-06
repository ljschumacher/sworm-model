% make movies from simulation results...

clear
close all

% general model parameters for all test - unless set otherwise
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
rc = 0.035;

revRatesClusterEdge = fliplr([0, 0.1, 0.2, 0.4, 0.8]);
speeds = [0.33];
slowspeeds = [0.33, 0.2, 0.1, 0.05];
attractionStrengths = [0];
paramCombis = combvec(revRatesClusterEdge,speeds,slowspeeds,attractionStrengths);
nParamCombis = size(paramCombis,2);
for paramCtr = 1:nParamCombis % can be parfor
    revRateClusterEdge = paramCombis(1,paramCtr);
    speed = paramCombis(2,paramCtr);
    slowspeed = paramCombis(3,paramCtr);
    attractionStrength = paramCombis(4,paramCtr);
    filename = ['woids_v0_' num2str(speed,'%1.0e') ...
        '_vs_' num2str(slowspeed,'%1.0e') ...
        '_epsLJ_' num2str(attractionStrength,'%1.0e')...
        '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e')];
    if exist(['results/woids/' filename '.mat'],'file')...
            &&~exist(['movies/woids/' filename '.mp4'],'file')
        out = load(['results/woids/' filename '.mat']);
        animateWoidTrajectories(out.xyarray,['movies/woids/' filename],L,rc);
    end
end
