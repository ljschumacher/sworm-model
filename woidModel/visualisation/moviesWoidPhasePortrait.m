% make movies from simulation results...

clear
close all

% general model parameters for all test - unless set otherwise
L = 7.5*[1, 1]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
rc = 0.035;
N = 40;
revRatesClusterEdge = 2%5:-1:0;%fliplr([0, 0.1, 0.2, 0.4, 0.8, 1.6]);
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125]);
slowspeeds = [0.018];
slowingMode = 'stochastic_bynode';
eps_LJ = 1e-3;
%f_hapt = 0.2;

k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = 0.6%fliplr(0:0.2:1);%fliplr([0 1./[8 4 2 1]]);
paramCombis = combvec(revRatesClusterEdge,speeds,slowspeeds,dkdN_dwell_values);
nParamCombis = size(paramCombis,2);
numRepeats = 1;
for repCtr =1:numRepeats
    for paramCtr = 1:nParamCombis % can be parfor but might impair movie quality
        revRateClusterEdge = paramCombis(1,paramCtr);
        speed = paramCombis(2,paramCtr);
        slowspeed = paramCombis(3,paramCtr);
        dkdN_dwell = paramCombis(4,paramCtr);
        filename = ['woids_N_' num2str(N) '_L_' num2str(L(1)) ...
            '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
            '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
            '_dkdN_' num2str(dkdN_dwell)...
            '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e')...
            '_LJsoft' num2str(eps_LJ) ...
            ...'_haptotaxis_' num2str(f_hapt) ...
            '_run' num2str(repCtr)];
        resultspath = '../results/woids/mapping/';
        moviepath = '../movies/woidMovies/mappingMovies/';
        if exist([resultspath filename '.mat'],'file')...
                &&~exist([moviepath filename '.mp4'],'file')
            out = load([resultspath filename '.mat']);
            animateWoidTrajectories(out.xyarray,[moviepath filename],L,rc);
            close(gcf)
        elseif ~exist([resultspath filename '.mat'],'file')
            disp(['no results for ' filename])
        end
    end
end