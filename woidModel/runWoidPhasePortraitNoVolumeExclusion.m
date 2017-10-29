% run simulations of simplified woid model with single node per woid
% for various speeds, attractions strengths, reversal probabilities...

% issues/todo:

clear
close all

% general model parameters for all test - unless set otherwise
N = 40; % N: number of objects
M = 49; % M: number of nodes in each object
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
numRepeats = 1;

T = 500; % T: simulation duration
rc = 0.035;
% saveevery = round(1/2/param.dT);
paramAll.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
% -- slow-down parameters --
paramAll.vs = 0;% vs: speed when slowed down (default v0/3)
paramAll.slowingNodes = [1:M];% slowingNodes: which nodes register contact (default head and tail)
paramAll.slowingMode = 'stochastic';
paramAll.k_dwell = 0.0036;
paramAll.k_undwell = 1.1;
% -- Lennard-Jones parameters --
paramAll.r_LJcutoff = 5*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
paramAll.sigma_LJ = 2*rc;  % particle size for Lennard-Jones force
% -- volume exclusion
paramAll.rc = 0;
% -- undulation parameters --
paramAll.angleNoise = 0;

revRatesClusterEdge = fliplr([0, 0.2, 0.4, 0.8, 1.6, 3.2]);
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.025, 0.0125, 0.005, 0.001]);
slowspeeds = [0.018];
attractionStrengths = [0];
dkdN_dwell_values = [0 1./[8 4 2 1 0.5]];
paramCombis = combvec(revRatesClusterEdge,speeds,slowspeeds,attractionStrengths,dkdN_dwell_values);
nParamCombis = size(paramCombis,2);
for repCtr = 1:numRepeats
for paramCtr = 1:nParamCombis
    param = paramAll;
    revRateClusterEdge = paramCombis(1,paramCtr);
    param.revRateClusterEdge = revRateClusterEdge;
    speed = paramCombis(2,paramCtr);
    param.v0 = speed;
    param.dT = min(1/2,rc/param.v0/16); % dT: time step, scales other parameters such as velocities and rates
    param.saveEvery = round(1/4/param.dT);
    param.vs = paramCombis(3,paramCtr);
    attractionStrength = paramCombis(4,paramCtr);
    param.dkdN_dwell = paramCombis(5,paramCtr);
    if attractionStrength>0
        param.r_LJcutoff = 5*rc;
    else
        param.r_LJcutoff = -1; % don't need to compute attraction if it's zero
    end
    param.eps_LJ = attractionStrength;
    filename = ['woids_N_' num2str(N) '_L_' num2str(L(1)) '_noVolExcl' ... 
        ... '_angleNoise'...
        '_v0_' num2str(param.v0,'%1.0e') '_vs_' num2str(param.vs,'%1.0e') ...
        '_' param.slowingMode 'SlowDown' '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
        '_dkdN_' num2str(param.dkdN_dwell)...num2str(param.num_nbr_max_per_node)...
        '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
        '_revRateClusterEdge_' num2str(param.revRateClusterEdge,'%1.0e') ...
        '_run' num2str(repCtr)];
    if ~exist(['results/woids/' filename '.mat'],'file')&&isempty(dir(['results/woids/' filename '_running_on_*.mat']))
        disp(['running ' filename])
        % make a dummy file to mark that this sim is running on this computer
        [~, hostname] = system('hostname -s'); hostname = strrep(hostname,newline,'');
        tmp_filename = ['results/woids/' filename '_running_on_' hostname '.mat'];
        save(tmp_filename,'N','M','L','param')
        rng(repCtr) % set random seed to be the same for each simulation
        xyarray = runWoids(T,N,M,L,param);
        xyarray = single(xyarray); % save space by using single precision
        saveResults(['results/woids/' filename '.mat'],struct('xyarray',xyarray,'T',T,'N',N,'M',M,'L',L,'param',param))
        delete(tmp_filename)
    end
end
end
