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

reversalMode = 'density';
drdN_rev_values = 0:0.2:1;

speed = [0.33];
% slowspeeds = fliplr([0.33, 0.05, 0.025, 0.0125]);
slowspeed = 0.018;
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = 0:0.2:1;
dkdN_roam_values = 0:0.2:2;
% angleNoise = 1;
k_theta = 2;
% f_hapt = 0.5;

secondVariables = dkdN_dwell_values;
ndrevVals = numel(drdN_rev_values);
ndwellVals = numel(dkdN_dwell_values);
nroamVals = numel(dkdN_roam_values);
aspectRatio = ndrevVals/(ndwellVals);

Rgyr = NaN(ndrevVals,ndwellVals,nroamVals); % radius of gyration
ri = 3*0.035;
%%
phasePortraitFig = figure;
for revRateCtr = 1:ndrevVals
    drdN_rev = drdN_rev_values(revRateCtr);
    for ddwellCtr = 1:ndwellVals
        dkdN_dwell = dkdN_dwell_values(ddwellCtr);
        for droamCtr = 1:nroamVals
            dkdN_undwell = dkdN_roam_values(droamCtr);
            filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
                ...'_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
                '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
                '_dkdN_' num2str(dkdN_dwell) '_' num2str(dkdN_undwell)... %% also check the case of indem dwell/roam
                '_rev' reversalMode '_drdN_' num2str(drdN_rev) ...
                ...'_haptotaxis_' num2str(f_hapt) ...
                '_clusteredStart' ...
                '_run1.mat'];
            filepath = '../results/woidlinos/';
            if exist([filepath filename],'file')
                load([filepath filename])
                time2plot = size(xyarray,4);
                positions2plot = xyarray(:,:,:,time2plot);
                % compute radius of gyration (of worm heads)
                Rgyr(revRateCtr,ddwellCtr,droamCtr) = sqrt(sum(var(positions2plot(:,1,:))));
            else
                warning([filename ' does not exist'])
            end
        end
    end
end
% make diagonal slice to get coordinates for subsequent plotting
dslice = surf(dkdN_dwell_values,drdN_rev_values,dkdN_dwell_values.*ones(size(drdN_rev_values))');
xd = dslice.XData;
yd = dslice.YData;
zd = dslice.ZData;
delete(dslice)

h = slice(dkdN_dwell_values,drdN_rev_values,dkdN_roam_values,log2(Rgyr),[0],[0.2 0],[ 0]);
hold on
h2 = slice(dkdN_dwell_values,drdN_rev_values,dkdN_roam_values,log2(Rgyr),...
    repmat(xd,ndrevVals,1),repmat(yd',1,ndwellVals),zd);
h = [h; h2];
nSlices = length(h);
for slcCtr = 1:nSlices
    h(slcCtr).FaceColor = 'interp';
    h(slcCtr).EdgeColor = 'none';
end
load ~/Dropbox/Utilities/colormaps_ascii/increasing_warm/cmap_RdOrYl.txt
colormap(flipud(cmap_RdOrYl(1:51:end,:)))
view(150,50)
ylabel('dr_{rev}/d\rho')
xlabel('dk_{dwell}/d\rho')
zlabel('dk_{roam}/d\rho')
% plot isosurface?
% isosurface(dkdN_dwell_values,drdN_rev_values,dkdN_roam_values,log2(Rgyr),1)
% add colorbar
hc = colorbar;
hc.Label.String = 'log_{2}R_{gyr}';
caxis([0 floor(log2(max(Rgyr(:))))])
%% export figure
phasePortraitFig.PaperUnits = 'centimeters';
filename = ['../figures/woidlinos/woidlinoPhaseDiagram_clusterStability_3D'...
    '_N_' num2str(N) '_M_' num2str(M) '_L_' num2str(L(1)) '_noVolExcl' ...
    ...'_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
    '_speed_' num2str(speed,'%1.0e') ...
    '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
    '_revDensity.eps'];
exportfig(phasePortraitFig,filename, exportOptions)
system(['epstopdf ' filename]);
system(['rm ' filename]);