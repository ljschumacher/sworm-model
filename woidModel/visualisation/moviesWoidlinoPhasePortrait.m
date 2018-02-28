% make movies from simulation results...

clear
close all

% general model parameters for all test - unless set otherwise
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
% rc = 0.035;
N = 40;
M = 18;
drdN_rev_values = fliplr(0.4:0.2:1);
speed = [0.33];
slowspeed = [0.018];
attractionStrengths = [0];
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = 0:0.2:1;
paramCombis = combvec(drdN_rev_values,dkdN_dwell_values);
nParamCombis = size(paramCombis,2);
% angleNoise = 0.02;

numRepeats = 1;
for repCtr =1:numRepeats
for paramCtr = 1:nParamCombis % can be parfor but might impair movie quality
    drdN_rev = paramCombis(1,paramCtr);
    dkdN_dwell = paramCombis(2,paramCtr);
    filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
        ...'_noVolExcl' '_angleNoise_' num2str(angleNoise) ...
        '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
        '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
        '_dkdN_' num2str(dkdN_dwell) ...
        ...'_epsLJ_' num2str(attractionStrength,'%1.0e') ...
        '_revdensity_drdN_' num2str(drdN_rev) ...
        '_run' num2str(repCtr)];
    resultspath = '../results/woidlinos/mapping/';
    moviepath = '../movies/woidlinoMovies/mappingMovies/';
    if exist([resultspath filename '.mat'],'file')...
            &&~exist([moviepath filename '.mp4'],'file')
        out = load([resultspath filename '.mat']);
        animateWoidTrajectories(out.xyarray,[moviepath filename],L);%,out.param.rc);
    elseif ~exist([resultspath filename '.mat'],'file')
        disp(['no results for ' filename])
    end
end
end