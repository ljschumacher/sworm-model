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
L = [2.4 2.4];
N = 2;
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
haptotaxisMode = 'weighted_additive';% -- speed and time-step --
v0 = 0.33; % npr1 0.33; N2 0.14

drdN_rev_values = linspace(0,1,10);
dkdN_dwell_values = linspace(0,1,10);
dkdN_undwell_values = linspace(0,2,10);
f_hapt_values = linspace(0,0.02,10);

filepath = '../results/woidlinos/pairedStart/';
filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta) ...
    '_v0_' num2str(v0,'%1.0e') '_vs_' num2str(vs,'%1.0e') ...
    '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
    '_rev' reversalMode ...
    '_haptotaxis_' haptotaxisMode ...
    '_pairedStart' ...
    '_minDistances.mat'];

load([filepath filename])

ndrevVals = numel(drdN_rev_values);
ndwellVals = numel(dkdN_dwell_values);
nundwellVals = numel(dkdN_undwell_values);

%% plot data
phasePortraitFig = figure;

% make diagonal slice to get coordinates for subsequent plotting
dslice = surf(dkdN_dwell_values,drdN_rev_values,...
    (-1.1+(fliplr(dkdN_undwell_values)'*(2+fliplr(dkdN_dwell_values))+0.5*fliplr(dkdN_dwell_values))));
xd = dslice.XData;
yd = dslice.YData;
zd = dslice.ZData;
delete(dslice)

% take median over replicate simulations
mminPdist = nanmedian(minPdist,5);

h = slice(dkdN_dwell_values,drdN_rev_values,dkdN_undwell_values,log2(mminPdist),[1],[],[ 0]);
hold on
h2 = slice(dkdN_dwell_values,drdN_rev_values,dkdN_undwell_values,log2(mminPdist),...
    repmat(xd,ndrevVals,1),repmat(yd',1,ndwellVals),zd);
h = [h; h2];
nSlices = length(h);
for slcCtr = 1:nSlices
    h(slcCtr).FaceColor = 'interp';
    h(slcCtr).EdgeColor = 'none';
end

% plot isosurface?
hs = patch(isosurface(dkdN_dwell_values,drdN_rev_values,dkdN_undwell_values,...
    log2(mminPdist),0));
hs.FaceAlpha = 0.5;

load ~/Dropbox/Utilities/colormaps_ascii/increasing_warm/cmap_RdOrYl.txt
colormap((cmap_RdOrYl(1:32:end,:)))
view(-45,30)
ylabel('dr_{rev}/d\rho')
xlabel('dk_{dwell}/d\rho')
zlabel('dk_{roam}/d\rho')
zlim([0 2])
%  add colorbar
hc = colorbar;
hc.Label.String = 'log_2<{R_{min}}>';
caxis([-4 floor(log2(max(mminPdist(:))))])
%% export figure
phasePortraitFig.PaperUnits = 'centimeters';
filename = ['../figures/woidlinos/woidlinoPhaseDiagram_pairStability_3D'...
    '_N_' num2str(N) '_M_' num2str(M) '_L_' num2str(L(1)) '_noVolExcl' ...
    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
    '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
    '_rev' reversalMode '_haptotaxis_' haptotaxisMode ...
    '_revDensity.eps'];
exportfig(phasePortraitFig,filename, exportOptions)
system(['epstopdf ' filename]);
system(['rm ' filename]);