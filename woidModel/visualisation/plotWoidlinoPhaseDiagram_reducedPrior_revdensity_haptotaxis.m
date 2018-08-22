% plot woidlino phase portrait
% shows the end-points of of woidlet simulations for a 2D parameter sweep
close all
clear

addpath('../')
addpath('../analysis/')

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',17,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1,...
    'Renderer','painters');

M = 18;
N = 40;
plotColor = [0.25, 0.25, 0.25];

% -- slow-down parameters --
vs = 0.018;% vs: speed when slowed down (default v0/3)
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
% -- reversal parameters --
reversalMode = 'density';
% -- Lennard-Jones parameters --
r_LJcutoff = -1;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
sigma_LJ = 0;  % particle size for Lennard-Jones force
eps_LJ = 0;
% -- undulation parameters --
k_theta = 0;
theta_0 = 0;
omega_m = 0;
deltaPhase = 0;
angleNoise = 0.05;
% -- haptotaxis --
Rif = 1.2/0.035;
haptotaxisMode = 'weighted_additive';
% -- speed and time-step --
v0 = 0.33; % npr1 0.33; N2 0.14

drdN_rev_values = linspace(0,1,10);
dkdN_dwell_values = linspace(0,1,10);
dkdN_undwell_values = linspace(0,2,10);
f_hapt_values = linspace(0,0.12,10);

ndrevVals = numel(drdN_rev_values);
ndwellVals = numel(dkdN_dwell_values);
nundwellVals = numel(dkdN_undwell_values);
nf_haptVals = numel(f_hapt_values);

supportLimits = [minmax(drdN_rev_values); minmax(dkdN_dwell_values);...
    minmax(dkdN_undwell_values); minmax(f_hapt_values)]';
paramCombis = combvec(drdN_rev_values,dkdN_dwell_values,dkdN_undwell_values,f_hapt_values)';

Rgyr = NaN(ndrevVals,ndwellVals,nundwellVals,nf_haptVals); % radius of gyration

pStabThresh = 1; % distance that pairs need to be apart
cStabThresh = 3; % min radius of gyration for clusters to be considered stable

%% load pair stability data
filepath = '../results/woidlinos/pairedStart/';
Lp = [2.4 2.4];
filename = ['wlM' num2str(M) '_N_2_L_' num2str(Lp(1)) ...
    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta) ...
    '_v0_' num2str(v0,'%1.0e') '_vs_' num2str(vs,'%1.0e') ...
    '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
    '_rev' reversalMode ...
    '_haptotaxis_' haptotaxisMode ...
    '_pairedStart' ...
    '_minDistances.mat'];

load([filepath filename])
% take median over replicate simulations
mminPdist = nanmedian(minPdist,5);

%% assemble cluster stability data
L = 1.8;
for revRateCtr = 1:ndrevVals
    drdN_rev = drdN_rev_values(revRateCtr);
    for ddwellCtr = 1:ndwellVals
        dkdN_dwell = dkdN_dwell_values(ddwellCtr);
        for dundwellCtr = 1:nundwellVals
            dkdN_undwell = dkdN_undwell_values(dundwellCtr);
            for f_haptCtr = 1:nf_haptVals
                f_hapt = f_hapt_values(f_haptCtr);
                filename = ['wlM' num2str(M) '_N_40_L_' num2str(L(1)) ...
                    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta) ...
                    '_v0_' num2str(v0,'%1.0e') '_vs_' num2str(vs,'%1.0e') ...
                    '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
                    '_dkdN_' num2str(dkdN_dwell) '_' num2str(dkdN_undwell)...
                    '_rev' reversalMode '_drdN_' num2str(drdN_rev) ...
                    '_haptotaxis_' haptotaxisMode '_' num2str(f_hapt) ...
                    '_clusteredStart' ...
                    '_run1.mat'];
                filepath = '../results/woidlinos/clusteredStart/npr_1/';
                if mminPdist(revRateCtr,ddwellCtr,dundwellCtr,f_haptCtr)>=pStabThresh...
                        &&exist([filepath filename],'file')
                    load([filepath filename])
                    % compute radius of gyration (of worm heads)
                    Rgyr(revRateCtr,ddwellCtr,dundwellCtr,f_haptCtr) = ...
                        sqrt(sum(var(xyarray(:,1,:))));
                elseif mminPdist(revRateCtr,ddwellCtr,dundwellCtr,f_haptCtr)>=pStabThresh...
                        &&~exist([filepath filename],'file')
                    warning([filename ' does not exist'])
                end
            end
        end
    end
end

%% kernel density estimation
% find points that are pair-unstable but cluster-stable
sampleIndcs = find(mminPdist(:)>=pStabThresh&Rgyr(:)<=cStabThresh);
nSamples = numel(sampleIndcs)

ndims = size(paramCombis,2);
bw = std(paramCombis(sampleIndcs,:)).*(4 + (ndims + 2)./nSamples).^(1./(ndims + 4)); %Silverman's rule of thumb for the bandwidth
% impose minimum bandwidth limit - set minimum sigma to be half grid size
if any(bw==0)
disp('adjusting zero bandwidths')
bw = max(bw,[mean(diff(drdN_rev_values)),...
mean(diff(dkdN_dwell_values)),...
mean(diff(dkdN_undwell_values)),...
mean(diff(f_hapt_values))]);
end

% % replicate multivariate kernel density estimation using a gaussian mixture model
% prior = mvksdensity(paramCombis(sampleIndcs,:),paramCombis,'Bandwidth',bw);
prior_npr1 = gmdistribution(paramCombis(sampleIndcs,:),bw.^2);

%% save prior for later sampling
save(['../priors4D_M_' num2str(M) '_noVolExcl' ...
    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
    '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
    '_rev' reversalMode '_haptotaxis_' haptotaxisMode ...
    '.mat'],'bw','cStabThresh','pStabThresh','supportLimits','prior_npr1')
%% plot marginals of reduced prior
priorFig = figure;
load('~/Dropbox/Utilities/colormaps_ascii/increasing_cool/cmap_Blues.txt')
% generate samples from prior for plotting
plotData = random(prior_npr1,1e6);
% truncate boundaries
plotData = plotData(all(plotData>0,2),:);
for dimCtr = 1:ndims
    overLogIndcs = plotData(:,dimCtr)>=supportLimits(2,dimCtr);
    plotData(overLogIndcs,:) = [];
end
[~, AX] = hplotmatrix(plotData,[],ones(size(plotData,1),1),supportLimits);
AX(1,1).YLabel.String = 'P';
AX(1,1).YTick = [];
AX(2,1).YLabel.String = 'dk_{dwell}/d\rho';
AX(2,1).YTick = AX(2,1).YTick(2:end); % remove overlapping tick-labels
AX(3,1).YLabel.String = 'dk_{roam}/d\rho';
AX(3,1).YTick = AX(3,1).YTick(2:end); % remove overlapping tick-labels
AX(4,1).YLabel.String = 'f_{taxis}';
AX(4,1).XLabel.String = 'dr_{rev}/d\rho';
AX(4,2).XLabel.String = 'dk_{dwell}/d\rho';
AX(4,3).XLabel.String = 'dk_{roam}/d\rho';
AX(4,4).XLabel.String = 'f_{taxis}';
colormap(flipud(cmap_Blues(1:24:end,:)))

%% export figure
priorFig.PaperUnits = 'centimeters';
filename = ['../figures/woidlinos/woidlinoPhaseDiagram_reducedPrior_4D'...
    '_M_' num2str(M) '_noVolExcl' ...
    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
    '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
    '_rev' reversalMode '_haptotaxis_' haptotaxisMode ...
    '.eps'];
exportfig(priorFig,filename, exportOptions)
system(['epstopdf ' filename]);
system(['rm ' filename]);