% plot woidlet phase portrait
% shows the end-points of of woidlet simulations for a 2D parameter sweep
close all
clear

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',17,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',6,...
    'LineWidth',1,...
    'Renderer','opengl');

radius = 0.035;
plotColor = [0.25, 0.25, 0.25];
Nval = 40;
Lval = 7.5;
revRatesClusterEdge = [0, 0.1, 0.2, 0.4, 0.8];
speeds = [0.33];
slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025]);
% slowspeeds = fliplr([0.33, 0.2, 0.1, 0.05]);
attractionStrength = [0];
for speed = speeds
    phasePortraitFig = figure;
    plotCtr = 1;
    phasePortraitFig.Name = ['speed = ' num2str(speed)];
    for slowspeed = slowspeeds
        for revRateClusterEdge = revRatesClusterEdge
            filename = ['../results/woids/woids_N_' num2str(Nval) '_L_' num2str(Lval) ...
                ...%'_noUndulations'...
                '_v0_' num2str(speed,'%1.0e') ...
                '_vs_' num2str(slowspeed,'%1.0e') '_gradualSlowDown' ...
                '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e') '.mat'];
            if exist(filename,'file')
                load(filename)
                time2plot = round(size(xyarray,4)*(0.9 + 0.1*rand()));
                positions2plot = xyarray(:,:,:,time2plot);
                subplot(length(slowspeeds),length(revRatesClusterEdge),plotCtr)
                ax = plotWoidTrajectoriesSingleFrame(positions2plot,L,radius,plotColor);
                title(['r=' num2str(revRateClusterEdge) ', v_s =' num2str(slowspeed,'%1.0e')],...
                    'FontWeight','normal')
                ax.Position = ax.Position.*[1 1 1.2 1.2] - [0.0 0.0 0 0]; % stretch panel
                ax.DataAspectRatio = [1 1 1];
                ax.Box = 'on';
            end
            plotCtr = plotCtr + 1;
        end
    end
    %% export figure
    phasePortraitFig.PaperUnits = 'centimeters';
    filename = ['../figures/woids/woidPhasePortrait_N_' num2str(Nval) '_L_' num2str(Lval) ...
        ...%'_noUndulations'...
        '_speed_' num2str(speed,'%1.0e') '_slowing' '_gradual' ...
        '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end