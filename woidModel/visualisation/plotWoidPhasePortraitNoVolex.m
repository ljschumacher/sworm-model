% plot woid phase portrait
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
    'Renderer','opengl');

radius = 0.035;
plotColor = [0.25, 0.25, 0.25];
N = 40;
M = 18;
Lval = 7.5;
revRatesClusterEdge = 0:5;
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125, 0.005]);
% slowspeeds = fliplr([0.33, 0.025, 0.0125, 0.005, 0.001]);
slowspeeds = [0.018];
slowingMode = 'stochastic_bynode';
% eps_LJ = 1e-3;
% f_hapt = 0.1;

k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = fliplr(0:0.2:1);
secondVariables = dkdN_dwell_values;
nrevRates = numel(revRatesClusterEdge);
ndwellVals = numel(dkdN_dwell_values);
aspectRatio = nrevRates/(ndwellVals + 1/3);
numRepeats = 1
% highlight panels
select_panels = [Inf, Inf]%[15, 13];
select_colors = lines(2);
for repCtr=1:numRepeats
    for speed = speeds
        phasePortraitFig = figure;
        plotCtr = 1;
        phasePortraitFig.Name = ['speed = ' num2str(speed)];
        for slowspeed = slowspeeds
            for dkdN_dwell = dkdN_dwell_values
                for revRateClusterEdge = revRatesClusterEdge
                    filename = ['woidsM' num2str(M) '_N_' num2str(N) '_L_' num2str(Lval) '_noVolExcl' ... 
                        '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                        '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
                        '_dkdN_' num2str(dkdN_dwell)...
                        '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e')...
                        ...'_LJsoft' num2str(eps_LJ) ...
                        ...'_haptotaxis_' num2str(f_hapt) ...
                        '_run' num2str(repCtr) '.mat'];
                    filepath = '../results/woids/mapping/';
                    if exist([filepath filename],'file')
                        load([filepath filename])
                        time2plot = round(size(xyarray,4)*(0.9 + 0.1*rand()));
                        positions2plot = xyarray(:,:,:,time2plot);
                        subplot(length(secondVariables),length(revRatesClusterEdge),plotCtr)
                        % highlight panels
                        if plotCtr==select_panels(1)
                            thisColor = select_colors(1,:);
                        elseif plotCtr==select_panels(2)
                            thisColor = select_colors(2,:);
                        else
                            thisColor = plotColor;
                        end
                        ax = plotWoidTrajectoriesSingleFrame(positions2plot,L,radius,thisColor,true);
                        ax.Position = ax.Position.*[1 1 1.25 1.25] - [0.0 0.0 0 0]; % stretch panel
                        ax.DataAspectRatio = [1 1 1];
                        ax.Box = 'on';
                    else
                        disp(['No file ' filename])
                    end
                    plotCtr = plotCtr + 1;
                end
            end
        end
        % make overall axes
        ax = axes('Color','none');
        ax.XTick = linspace(1/nrevRates/2,1-1/nrevRates/2,nrevRates);
        ax.XTickLabel = num2str(revRatesClusterEdge');
        ax.YTick = linspace(1/ndwellVals/2,1-1/ndwellVals/2,ndwellVals);
        ax.YTickLabel = num2str(fliplr(dkdN_dwell_values)');
        ax.TickDir = 'out';
        ax.Position = ax.Position.*[1 1 1 1];
        xlabel('cluster-edge reversal rate (1/s)')
        ylabel('dk/d\rho')
        %% export figure
        phasePortraitFig.Position(3) = phasePortraitFig.Position(4)*aspectRatio; % resize figure
        phasePortraitFig.PaperUnits = 'centimeters';
        filename = ['../figures/woids/woidPhasePortrait_mapping_' num2str(repCtr) '_N_' num2str(N) '_L_' num2str(Lval) ...
            '_noVolExcl'...'_noUndulations' ...'_angleNoise'
            '_speed_' num2str(speed,'%1.0e') '_slowing' '_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
            ...'_LJsoft' num2str(eps_LJ) ...
            ...'_haptotaxis_' num2str(f_hapt) ...
            '.eps'];
        exportfig(phasePortraitFig,filename, exportOptions)
        system(['epstopdf ' filename]);
        system(['rm ' filename]);
    end
end