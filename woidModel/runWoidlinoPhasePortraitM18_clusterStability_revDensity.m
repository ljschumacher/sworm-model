% run simulations of simplified woid model with single node per woid
% for various speeds, attractions strengths, reversal probabilities...
% start with clustered initial conditions, and free boundary, to assess
% parameter bounds on cluster stability

% issues/todo:

clear
close all

% general model parameters for all test - unless set otherwise
N = 40; % N: number of objects
M = 18; % M: number of nodes in each object
L = 1.8; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
numRepeats = 1;

T = 300;
rc0 = 0.035; % rc: core repulsion radius (default 0.035 mm)
param.rc = 0;
param.ri = 3*rc0;
param.bc = 'free'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
param.segmentLength = 1.13/(M - 1);
% -- slow-down parameters --
param.vs = 0.018;% vs: speed when slowed down (default v0/3)
param.slowingNodes = 1:M;% slowingNodes: which nodes register contact (default head and tail)
param.slowingMode = 'stochastic_bynode';
param.k_dwell = 0.0036;
param.k_undwell = 1.1;
% -- reversal parameters --
param.reversalMode = 'density';
param.revRateClusterEdge = 0;
% -- Lennard-Jones parameters --
param.r_LJcutoff = -1;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
param.sigma_LJ = 0;  % particle size for Lennard-Jones force
param.eps_LJ = 0;
% -- undulation parameters --
param.k_theta = 0;
param.theta_0 = 0;
param.omega_m = 0;
param.deltaPhase = 0;
param.angleNoise = 0.05;
% -- haptotaxis --
param.Rif = 1.2/0.035;
param.haptotaxisMode = 'weighted_additive';
% -- speed and time-step --
param.v0 = 0.33; % npr1 0.33; N2 0.14
param.dT = min(1/2,rc0/param.v0/8); % dT: time step, scales other parameters such as velocities and rates
param.saveEvery = round(1/param.dT);

drdN_rev_values = linspace(0,1,10);
dkdN_dwell_values = linspace(0,1,10);
dkdN_undwell_values = linspace(0,2,10);
f_hapt_values = linspace(0,0.02,10);

ndrevVals = numel(drdN_rev_values);
ndwellVals = numel(dkdN_dwell_values);
nundwellVals = numel(dkdN_undwell_values);
nf_haptVals = numel(f_hapt_values);

% load pair stability data
filepath = 'results/woidlinos/pairedStart/';
filename = ['wlM' num2str(M) '_N_2_L_2.4' ...
    '_angleNoise_' num2str(param.angleNoise) '_k_theta_' num2str(param.k_theta) ...
    '_v0_' num2str(param.v0,'%1.0e') '_vs_' num2str(param.vs,'%1.0e') ...
    '_' param.slowingMode 'SlowDown' '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
    '_rev' param.reversalMode ...
    '_haptotaxis_' param.haptotaxisMode ...
    '_pairedStart' ...
    '_minDistances.mat'];
load([filepath filename],'minPdist')
% take median over replicate simulations
minPdist = nanmedian(minPdist,5);

pStabThresh = 1; % distance that pairs need to be apart

for repCtr = 1:numRepeats
    for revRateCtr = 1:ndrevVals
        param.drdN_rev = drdN_rev_values(revRateCtr);
        for ddwellCtr = 1:ndwellVals
            param.dkdN_dwell = dkdN_dwell_values(ddwellCtr);
            for dundwellCtr = 1:nundwellVals
                param.dkdN_undwell = dkdN_undwell_values(dundwellCtr);
                for f_haptCtr = 1:nf_haptVals
                    param.f_hapt = f_hapt_values(f_haptCtr);
                    % only run this simulation if the parameter combination
                    % does not result in stable pairs
                    if minPdist(revRateCtr,ddwellCtr,dundwellCtr,f_haptCtr)...
                            >=pStabThresh
                        filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
                            '_angleNoise_' num2str(param.angleNoise) '_k_theta_' num2str(param.k_theta) ...
                            '_v0_' num2str(param.v0,'%1.0e') '_vs_' num2str(param.vs,'%1.0e') ...
                            '_' param.slowingMode 'SlowDown' '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell) ...
                            '_dkdN_' num2str(param.dkdN_dwell) '_' num2str(param.dkdN_undwell)...
                            '_rev' param.reversalMode '_drdN_' num2str(param.drdN_rev) ...
                            '_haptotaxis_' param.haptotaxisMode '_' num2str(param.f_hapt) ...
                            '_clusteredStart' ...
                            '_run' num2str(repCtr)];
                        filepath = 'results/woidlinos/clusteredStart/';
                        if ~exist([filepath filename '.mat'],'file')...
                                &&isempty(dir([filepath filename '_running_on_*.mat']))
                            disp(['running ' filename])
                            % make a dummy file to mark that this sim is running on this computer
                            [~, hostname] = system('hostname -s'); hostname = strrep(hostname,newline,'');
                            tmp_filename = [filepath filename '_running_on_' hostname '.mat'];
                            save(tmp_filename,'N','M','L','param')
                            rng(repCtr) % set random seed to be the same for each simulation
                            [xyarray, currentState] = runWoids(T,N,M,L,param);
                            xyarray = single(xyarray); % save space by using single precision
                            xyarray = xyarray(:,:,:,end); % save end-point only
                            save([filepath filename '.mat'],'xyarray','T','N','M','L','param','currentState')
                            delete(tmp_filename)
                        end
                    else
                        disp(['params result in stable pair (median min separation ' ...
                            num2str(minPdist(revRateCtr,ddwellCtr,dundwellCtr,f_haptCtr)) '), discontinuing simulation'])
                    end
                end
            end
        end
    end
end
