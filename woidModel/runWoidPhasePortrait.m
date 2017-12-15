function [] = runWoidPhasePortrait(N,L)
% run simulations of simplified woid model with single node per woid
% for various speeds, attractions strengths, reversal probabilities...

% issues/todo:

% general model parameters for all test - unless set otherwise
% N = 60; % N: number of objects
M = 49; % M: number of nodes in each object
% L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
if numel(L)==1
    L = [L, L];
end
T = 1000; % T: simulation duration
numRepeats = 1;
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
paramAll.r_LJcutoff = 4*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
paramAll.sigma_LJ = 2*rc;  % particle size for Lennard-Jones force

paramAll.LJmode = 'soft';
paramAll.r_LJcutoff = 2*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
paramAll.rc = 0;

revRatesClusterEdge = 6.4%fliplr([0, 0.4, 0.8, 1.6, 3.2, 6.4]);
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125, 0.005]);
slowspeeds = [0.018];
attractionStrengths = [1e-2];
% num_nbr_max_per_nodes = [3 4];
dkdN_dwell_values = [0.125]%[0 1./[8 4 2 1 0.5]];
paramCombis = combvec(revRatesClusterEdge,speeds,slowspeeds,attractionStrengths,dkdN_dwell_values);
nParamCombis = size(paramCombis,2);
for paramCtr = 1:nParamCombis
    param = paramAll;
    revRateClusterEdge = paramCombis(1,paramCtr);
    param.revRateClusterEdge = revRateClusterEdge;
    speed = paramCombis(2,paramCtr);
    param.omega_m = 2*pi*0.6/0.33*speed;
    param.v0 = speed;
    param.dT = min(1/2,rc/param.v0/16); % dT: time step, scales other parameters such as velocities and rates
    param.saveEvery = round(1/4/param.dT);
    param.vs = paramCombis(3,paramCtr);
    attractionStrength = paramCombis(4,paramCtr);
    param.dkdN_dwell = paramCombis(5,paramCtr);
    if attractionStrength==0
        param.r_LJcutoff = -1; % don't need to compute attraction if it's zero
    end
    param.eps_LJ = attractionStrength;
    for repCtr = 1:numRepeats
    filename = ['woids_N_' num2str(N) '_L_' num2str(L(1)) ...
        '_v0_' num2str(param.v0,'%1.0e') '_vs_' num2str(param.vs,'%1.0e') ...
        '_' param.slowingMode 'SlowDown' '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
        '_dkdN_' num2str(param.dkdN_dwell)...
        '_epsLJ_' num2str(attractionStrength,'%1.0e') '_' param.LJmode ...
        '_noContactForce' ...
        '_revRateClusterEdge_' num2str(param.revRateClusterEdge,'%1.0e')...
        '_run' num2str(repCtr)];
    if ~exist(['results/woids/' filename '.mat'],'file')&&isempty(dir(['results/woids/' filename '_running_on_*.mat']))
        disp(['running ' filename])
        % make a dummy file to mark that this sim is running on this computer
        [~, hostname] = system('hostname -s'); hostname = strrep(hostname,newline,'');
        tmp_filename = ['results/woids/' filename '_running_on_' hostname '.mat'];
        save(tmp_filename,'N','M','L','param')
        rng(1) % set random seed to be the same for each simulation
        [xyarray, currentState] = runWoids(T,N,M,L,param);
        xyarray = single(xyarray); % save space by using single precision
        saveResults(['results/woids/' filename '.mat'],...
            struct('xyarray',single(xyarray),'T',T,'N',N,'M',M,'L',L,'param',param,'currentState',currentState))
        delete(tmp_filename)
    end
    end
end
