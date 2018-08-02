% plot woidlet phase portrait
% shows the end-points of of woidlet simulations for a 2D parameter sweep
close all
clear

addpath('../')

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',17,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1,...
    'Renderer','opengl');

radius = 0.035;
plotColor = [0.25, 0.25, 0.25];
Nval = 40;
Lval = 7.5;
revRatesClusterEdge = [0, 0.2, 0.4, 0.8, 1.6];
speeds = [0.33];
slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125]);
attractionStrength = [0];
slowingMode = 'gradual';
secondVariables = slowspeeds;
nrevRates = numel(revRatesClusterEdge);
nslowpeeds = numel(slowspeeds);
aspectRatio = nrevRates/nslowpeeds;
for speed = speeds
    phasePortraitFig = figure;
    plotCtr = 1;
    phasePortraitFig.Name = ['speed = ' num2str(speed)];
    for slowspeed = slowspeeds
            for revRateClusterEdge = revRatesClusterEdge
                filename = ['../results/woids/woids_N_' num2str(Nval) '_L_' num2str(Lval) ...
                    '_noUndulations'...'_noVolExcl' ...'_angleNoise'
                    '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e')...
                    '_' slowingMode 'SlowDown' ...
                    '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                    '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e')...
                    '.mat'];
                if exist(filename,'file')
                    load(filename)
                    time2plot = round(size(xyarray,4));%*(0.9 + 0.1*rand()));
                    positions2plot = xyarray(:,:,:,time2plot);
                    subplot(length(secondVariables),length(revRatesClusterEdge),plotCtr)
                    ax = plotWoidTrajectoriesSingleFrame(positions2plot,L,radius,plotColor,true);
%                     title(['r=' num2str(revRateClusterEdge) ', dk/dN =' num2str(dkdN_dwell,'%1.0e')],...
%                         'FontWeight','normal')
                    ax.Position = ax.Position.*[1 1 1.25 1.25] - [0.0 0.0 0 0]; % stretch panel
                    ax.DataAspectRatio = [1 1 1];
                    ax.Box = 'on';
                end
                plotCtr = plotCtr + 1;
            end
    end
    % make overall axes
    ax = axes('Color','none');
    ax.XTick = linspace(1/nrevRates/2,1-1/nrevRates/2,nrevRates);
    ax.XTickLabel = num2str(revRatesClusterEdge');
    ax.YTick = linspace(1/nslowpeeds/2,1-1/nslowpeeds/2,nslowpeeds);
    ax.YTickLabel = num2str(fliplr(slowspeeds)');
    ax.TickDir = 'out';
    xlabel('cluster-edge reversal rate (1/s)')
    ylabel('dk/d\rho')
    %% export figure
    phasePortraitFig.Position(3) = phasePortraitFig.Position(4)*aspectRatio; % resize figure
    phasePortraitFig.PaperUnits = 'centimeters';
    filename = ['../figures/woids/woidPhasePortrait_N_' num2str(Nval) '_L_' num2str(Lval) ...
        '_noUndulations'...'_noVolExcl' ...'_angleNoise'
        '_speed_' num2str(speed,'%1.0e') '_slowing' '_' slowingMode ...
        '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end