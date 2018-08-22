% plot woidlino phase portrait
% shows the end-points of of woidlet simulations for a 2D parameter sweep
close all
clear

addpath('../')

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',17,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1,...
    'Renderer','painters');

M = 18;
L = [7.5 7.5];
N = 40;
plotColor = [0.2, 0.2, 0.2];


T = 1500;
rc0 = 0.035; % rc: core repulsion radius (default 0.035 mm)
rc = 0;
ri = 3*rc0;
bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
segmentLength = 1.13/(M - 1);
% -- slow-down parameters --
vs = 0.018;% vs: speed when slowed down (default v0/3)
slowingNodes = 1:M;% slowingNodes: which nodes register contact (default head and tail)
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
% -- reversal parameters --
reversalMode = 'density';
revRateClusterEdge = 0;
% -- Lennard-Jones parameters --
r_LJcutoff = -1;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
sigma_LJ = 0;  % particle size for Lennard-Jones force
eps_LJ = 0;
% -- undulation parameters --
theta_0 = 0;
omega_m = 0;
deltaPhase = 0;
angleNoise = 0.05;
k_theta = 0;
% -- haptotaxis --
% f_hapt = 0.5;
% -- speed and time-step --
v0 = [0.33]; % npr1 0.33; N2 0.14
dT = min(1/2,rc0/v0/8); % dT: time step, scales other parameters such as velocities and rates
saveEvery = round(1/dT);

drdN_rev_values = linspace(0,1,5);
dkdN_dwell_values = fliplr(linspace(0,1,5));

secondVariables = dkdN_dwell_values;
ndrevRates = numel(drdN_rev_values);
ndwellVals = numel(dkdN_dwell_values);
aspectRatio = ndrevRates/(ndwellVals);

% highlight panels
select_panels = [7];
select_colors = lines(1);
for repCtr =1:1
    phasePortraitFig = figure;
    plotCtr = 1;
    for dkdN_dwell = dkdN_dwell_values
        for drdN_rev = drdN_rev_values
            filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
                '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta) ...
                '_v0_' num2str(v0,'%1.0e') '_vs_' num2str(vs,'%1.0e') ...
                '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
                '_dkdN_' num2str(dkdN_dwell) '_' num2str(dkdN_dwell)...
                '_rev' reversalMode '_drdN_' num2str(drdN_rev) ...
                ...'_haptotaxis_' num2str(f_hapt) ...
                '_run' num2str(repCtr) '.mat'];
            filepath = '../results/woidlinos/floppy/';
            if exist([filepath filename],'file')
                load([filepath filename])
                time2plot = round(size(xyarray,4)*(0.9 + 0.1*rand()));
                positions2plot = xyarray(:,:,:,time2plot);
                subplot(length(secondVariables),length(drdN_rev_values),plotCtr)
                % highlight panels
                %                         if plotCtr==select_panels(1)
                %                             thisColor = select_colors(1,:);
                %                         elseif plotCtr==select_panels(2)
                %                             thisColor = select_colors(2,:);
                if ismember(plotCtr,select_panels)
                    thisColor = select_colors;
                else
                    thisColor = plotColor;
                end
                radius = max(0.035,param.rc);
                ax = plotWoidTrajectoriesSingleFrame(positions2plot,L,radius,thisColor,true);
                ax.Position = ax.Position.*[1 1 1.25 1.25] - [0.0 0.0 0 0]; % stretch panel
                ax.DataAspectRatio = [1 1 1];
                ax.Box = 'on';
            else
                warning([filename ' does not exist'])
            end
            plotCtr = plotCtr + 1;
        end
    end
    % make overall axes
    ax = axes('Color','none');
    ax.XTick = linspace(1/ndrevRates/2,1-1/ndrevRates/2,ndrevRates);
    ax.XTickLabel = num2str(drdN_rev_values');
    ax.YTick = linspace(1/ndwellVals/2,1-1/ndwellVals/2,ndwellVals);
    ax.YTickLabel = num2str(fliplr(dkdN_dwell_values)');
    ax.TickDir = 'out';
    ax.Position = ax.Position.*[1 1 1 1];
    xlabel('dr_{rev}/d\rho')
    ylabel('dk/d\rho')
    %% export figure
    phasePortraitFig.Position(3) = phasePortraitFig.Position(4)*aspectRatio; % resize figure
    phasePortraitFig.PaperUnits = 'centimeters';
    filename = ['../figures/woidlinos/woidlinoPhasePortrait_mapping_' num2str(repCtr) '_N_' num2str(N) ...
        '_M_' num2str(M) '_L_' num2str(L(1)) '_noVolExcl' ...
        '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
        '_vo_' num2str(v0,'%1.0e') ...
        '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
        '_rev' reversalMode...'_haptotaxis_' num2str(f_hapt) ...
        '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end
