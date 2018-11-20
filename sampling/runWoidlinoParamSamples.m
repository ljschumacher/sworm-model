function [] = runWoidlinoParamSamples(sampleCtr)
% run simulations of simplified woid model with single node per woid
% for previously generated random parameter samples
addpath('../')

% general model parameters for all simulations - unless set otherwise
N = 40; % N: number of objects
M = 18; % M: number of nodes in each object
L = [7.5, 7.5]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
T = 2000;
rc = 0.035; % rc: core repulsion radius (default 0.035 mm)
param.rc = 0;
param.ri = 3*rc;
param.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
param.segmentLength = 1.13/(M - 1);
% -- slow-down parameters --
param.vs = 0.018; % npr1 0.018; N2 0.014
param.slowingNodes = 1:M;% slowingNodes: which nodes register contact (default head and tail)
param.slowingMode = 'stochastic_bynode';
param.k_dwell = 0.0036; % npr1 0.0036; N2 0.25
param.k_undwell = 1.1; % npr1 1.1; N2 0.45
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
param.v0 = [0.33]; % npr1 0.33; N2 0.14
param.dT = min(1/2,rc/param.v0/8); % dT: time step, scales other parameters such as velocities and rates
param.saveEvery = round(1/param.dT);

% load parameter samples
load(['paramSamples_nSamples100000_log_PRW_4D_wa_r2_npr1'...
    '.mat'],'paramSamples','supportLimits')
% set model parameters from generated samples
param.drdN_rev = paramSamples.drdN_rev(sampleCtr);
param.dkdN_dwell = paramSamples.dkdN_dwell(sampleCtr);
param.dkdN_undwell = paramSamples.dkdN_undwell(sampleCtr);
param.f_hapt = paramSamples.f_hapt(sampleCtr);

% filepath = '/exports/eddie/scratch/lschuma2/woidlinos/PRW_4D_r2/npr_1/';
filepath = '../results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r2/npr_1/';
filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
    '_v0_' num2str(param.v0) '_vs_' num2str(param.vs) ...
    '_angleNoise_' num2str(param.angleNoise) '_k_theta_' num2str(param.k_theta)...
    '_slow_' param.slowingMode '_dwell_' num2str(param.k_dwell) '_' num2str(param.k_undwell)...
    '_dkdN_' num2str(param.dkdN_dwell) '_' num2str(param.dkdN_undwell)...
    '_rev' param.reversalMode '_drdN_' num2str(param.drdN_rev) ...
    '_haptotaxis_' param.haptotaxisMode '_' num2str(param.f_hapt) ...
    '_sample_' num2str(sampleCtr)];
pStabThresh = 1;
cStabThresh = 3;
if ~exist([filepath filename '.mat'],'file')
    %% check for pair stability, which we don't want
    disp('Checking pair stability...')
    % set up paired initial conditions
    rng(5) % this happens to give a good pair of initial positions
    param.bc = 'free';
    [~, initialState] = runWoids(1,2,M,[2.4, 2.4],param);
    minPdist = zeros(10,1);
    for repCtr = 1:10
        rng(repCtr)
        [pairxyarray, ~] = runWoids(60,2,M,[2.4, 2.4],param,'resumeState',initialState);
        % compute minimum seperation
        worm1 = squeeze(pairxyarray(1,:,:,end));
        worm2 = squeeze(pairxyarray(2,:,:,end));
        X = worm1(:,1) - worm2(:,1)';
        Y = worm1(:,2) - worm2(:,2)';
        D = sqrt(X.^2 + Y.^2);
        minPdist(repCtr) = min(D(:));
    end
    if median(minPdist)<pStabThresh
        disp(['params result in stable pair (median min separation ' ...
            num2str(median(minPdist)) '), discontinuing simulation'])
        % save minPdist result?
    elseif median(minPdist)>=pStabThresh
        %% check for cluster stability, which we do/don't want (npr1/N2)
        disp('Checking cluster stability...')
        param.bc = 'free';
        rng(sampleCtr)
        [clustxyarray, ~] = runWoids(300,N,M,1.8,param);
        % compute radius of gyration (of worm heads)
        Rgyr = sqrt(sum(var(clustxyarray(:,1,:,end))));
        if Rgyr>cStabThresh
            disp(['params result in unstable cluster (Rgyr ' num2str(Rgyr) ...
                '), discontinuing simulation'])
            % save Rgyr result?
        elseif Rgyr<=cStabThresh
            %% run full-length simulation
            disp('Running full-length simulation...')
            param.bc = 'periodic';
            rng(sampleCtr) % set random seed to be DIFFERENT for each simulation
            [xyarray, currentState] = runWoids(T,N,M,L,param);
            xyarray = single(xyarray); % save space by using single precision
            save([filepath filename '.mat'],'xyarray','T','N','M','L','param','currentState')
        end
    end
end
end
