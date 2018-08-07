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
L = [3.6 3.6];
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
% f_hapt = 0.5;
% -- speed and time-step --
v0 = 0.33; % npr1 0.33; N2 0.14

drdN_rev_values = 0:0.1:1;
dkdN_dwell_values = 0:0.1:1;
dkdN_undwell_values = 0:0.2:2;

secondVariables = dkdN_dwell_values;
ndrevVals = numel(drdN_rev_values);
ndwellVals = numel(dkdN_dwell_values);
nundwellVals = numel(dkdN_undwell_values);
aspectRatio = ndrevVals/(ndwellVals);

Rgyr = NaN(ndrevVals,ndwellVals,nundwellVals); % radius of gyration
ri = 3*0.035;
%%
for revRateCtr = 1:ndrevVals
    drdN_rev = drdN_rev_values(revRateCtr);
    for ddwellCtr = 1:ndwellVals
        dkdN_dwell = dkdN_dwell_values(ddwellCtr);
        for dundwellCtr = 1:nundwellVals
            dkdN_undwell = dkdN_undwell_values(dundwellCtr);
            filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
                '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta) ...
                '_v0_' num2str(v0,'%1.0e') '_vs_' num2str(vs,'%1.0e') ...
                '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
                '_dkdN_' num2str(dkdN_dwell) '_' num2str(dkdN_undwell)...
                '_rev' reversalMode '_drdN_' num2str(drdN_rev) ...
                ...'_haptotaxis_' num2str(f_hapt) ...
                '_clusteredStart' ...
                '_run1.mat'];
            filepath = '../results/woidlinos/clusteredStart/';
            if exist([filepath filename],'file')
                load([filepath filename])
                time2plot = size(xyarray,4);
                positions2plot = xyarray(:,:,:,time2plot);
                % compute radius of gyration (of worm heads)
                Rgyr(revRateCtr,ddwellCtr,dundwellCtr) = sqrt(sum(var(positions2plot(:,1,:))));
            else
                warning([filename ' does not exist'])
            end
        end
    end
end
%% make diagonal slice to get coordinates for subsequent plotting
phasePortraitFig = figure;
dslice = surf(dkdN_dwell_values,drdN_rev_values,ones(size(drdN_rev_values))...
    .*fliplr(dkdN_undwell_values.^2)'./2);
xd = dslice.XData;
yd = dslice.YData;
zd = dslice.ZData;
delete(dslice)
% % make another slice
dslice3 = surf(dkdN_dwell_values,drdN_rev_values,(fliplr(dkdN_undwell_values.^2) ...
        .*fliplr(drdN_rev_values.^2)'));
xd3 = dslice3.XData;
yd3 = dslice3.YData;
zd3 = dslice3.ZData;
delete(dslice3)

h = slice(dkdN_dwell_values,drdN_rev_values,dkdN_undwell_values,log2(Rgyr),[],[0.1],[]);
hold on
h2 = slice(dkdN_dwell_values,drdN_rev_values,dkdN_undwell_values,log2(Rgyr),...
    repmat(xd,ndrevVals,1),repmat(yd',1,ndwellVals),zd-0.2);
h3 = slice(dkdN_dwell_values,drdN_rev_values,dkdN_undwell_values,log2(Rgyr),...
    repmat(xd3,ndrevVals,1),repmat(yd3',1,ndwellVals),zd3);
h = [h; h2; h3];
nSlices = length(h);
for slcCtr = 1:nSlices
    h(slcCtr).FaceColor = 'interp';
    h(slcCtr).EdgeColor = 'none';
end
load ~/Dropbox/Utilities/colormaps_ascii/increasing_warm/cmap_RdOrYl.txt
colormap(flipud(cmap_RdOrYl(1:51:end,:)))
view(120,30)
zlim([0 2])
ylabel('dr_{rev}/d\rho')
xlabel('dk_{dwell}/d\rho')
zlabel('dk_{roam}/d\rho')
% plot isosurface?
hs = patch(isosurface(dkdN_dwell_values,drdN_rev_values,dkdN_undwell_values,...
    log2(Rgyr),2));
hs.FaceAlpha = 0.5;
% add colorbar
hc = colorbar;
hc.Label.String = 'log_{2}R_{gyr}';
caxis([0 floor(log2(max(Rgyr(:))))])
%% export figure
phasePortraitFig.PaperUnits = 'centimeters';
filename = ['../figures/woidlinos/woidlinoPhaseDiagram_clusterStability_3D'...
    '_N_' num2str(N) '_M_' num2str(M) '_L_' num2str(L(1)) '_noVolExcl' ...
    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
    '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
    '_rev' param.reversalMode ...'_haptotaxis_' num2str(f_hapt) ...
    '_revDensity.eps'];
exportfig(phasePortraitFig,filename, exportOptions)
system(['epstopdf ' filename]);
system(['rm ' filename]);