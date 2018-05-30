function [] = runWoidPhasePortraitMidRangeRevDensity(paramCtr)
% run simulations of simplified woid model with single node per woid
% for various speeds, attractions strengths, reversal probabilities...

% issues/todo:

% general model parameters for all test - unless set otherwise
N = 40; % N: number of objects
M = 36; % M: number of nodes in each object
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
if numel(L)==1
    L = [L, L];
end
T = 1000; % T: simulation duration
rc = 0.035;
paramAll.ri = 7.5*rc;
% saveevery = round(1/2/param.dT);
paramAll.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
paramAll.segmentLength = 1.13/(M - 1);
paramAll.k_l = 80;
% -- slow-down parameters --
paramAll.vs = 0.018;% vs: speed when slowed down (default v0/3)
paramAll.slowingNodes = [1:M];% slowingNodes: which nodes register contact (default head and tail)
paramAll.slowingMode = 'stochastic_bynode';
paramAll.k_dwell = 0.0036;
paramAll.k_undwell = 1.1;
% -- reversal parameters --
paramAll.reversalMode = 'density';
paramAll.revRateClusterEdge = 0;
% -- Lennard-Jones parameters --
paramAll.r_LJcutoff = 3.75*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
paramAll.sigma_LJ = 2*rc;  % particle size for Lennard-Jones force
paramAll.eps_LJ = 0;%1e-3;
if paramAll.eps_LJ<=0
    paramAll.r_LJcutoff = -1; % don't need to compute attraction if it's zero
end
paramAll.LJmode = 'soft';
% % -- haptotaxis
paramAll.f_hapt = 0;%0.1;
% -- speed and time-step --
paramAll.v0 = [0.33]; % npr1 0.33; N2 0.14
paramAll.dT = min(1/2,rc/paramAll.v0/16); % dT: time step, scales other parameters such as velocities and rates
paramAll.saveEvery = round(1/paramAll.dT);

numRepeats = 1;
drdN_rev_values = 0:0.025:0.125;%0.05:0.25;
dkdN_dwell_values = 0:0.025:0.125;%:0.05:0.25;
dkdN_undwell_values = 0:0.025:0.25;%:0.05:0.5;

paramCombis = combvec(drdN_rev_values,dkdN_dwell_values,dkdN_undwell_values);
for repCtr = 1:numRepeats
%     for paramCtr = 1:nParamCombis
        param = paramAll;
        param.drdN_rev =  paramCombis(1,paramCtr);
        param.dkdN_dwell = paramCombis(2,paramCtr);
        param.dkdN_undwell = paramCombis(3,paramCtr);
        filename = ['woids_N_' num2str(N) '_L_' num2str(L(1)) ...
            '_v0_' num2str(param.v0,'%1.0e') '_vs_' num2str(param.vs,'%1.0e') ...
            '_' param.slowingMode 'SlowDown' '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
            '_dkdN_' num2str(param.dkdN_dwell) '_' num2str(param.dkdN_undwell)...
            '_rev' param.reversalMode '_drdN_' num2str(param.drdN_rev) ...
            ...'_LJ' param.LJmode num2str(param.eps_LJ) ...
            ...'_haptotaxis_' num2str(param.f_hapt) ...
            '_ri' num2str(param.ri) ...
            '_run' num2str(repCtr)];
%         filepath = 'results/woids/mapping/';
        filepath = '/work/lschumac/woids/';
        if ~exist([filepath filename '.mat'],'file')%...
% %                 &&isempty(dir([filepath filename '_running_on_*.mat']))
%             disp(['running ' filename])
%             % make a dummy file to mark that this sim is running on this computer
%             [~, hostname] = system('hostname -s'); hostname = strrep(hostname,newline,'');
%             tmp_filename = [filepath filename '_running_on_' hostname '.mat'];
%             save(tmp_filename,'N','M','L','param')
            rng(repCtr) % set random seed to be the same for each simulation
            [xyarray, currentState] = runWoids(T,N,M,L,param);
            xyarray = single(xyarray); % save space by using single precision
            save([filepath filename '.mat'],'xyarray','T','N','M','L','param','currentState')
%             delete(tmp_filename)
        end
%     end
end
end
