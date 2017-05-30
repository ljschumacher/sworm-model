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
    'LineWidth',1);

radius = 0.035;
plotColor = [0.25, 0.25, 0.25];

revRatesClusterEdge = [0, 0.1, 0.2, 0.4, 0.8];
speeds = [0.33, 0.14];
attractionStrengths = [0, 1e-5, 5e-5, 1e-4];
for speed = speeds
    phasePortraitFig = figure;
    plotCtr = 1;
    phasePortraitFig.Name = ['speed = ' num2str(speed)];
    for attractionStrength = attractionStrengths
        for revRateClusterEdge = revRatesClusterEdge
            filename = ['results/woids/woids_v0_' num2str(speed,'%1.0e') ...
                '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e') '.mat'];
            if exist(filename,'file')
                load(filename)
                positions2plot = xyarray(:,:,:,end);
                subplot(length(attractionStrengths),length(revRatesClusterEdge),plotCtr)
                ax = plotWoidTrajectoriesSingleFrame(positions2plot,L,radius);%,plotColor);
                title(['r=' num2str(revRateClusterEdge) ', \epsilon =' num2str(attractionStrength,'%1.0e')],...
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
    filename = ['figures/woids/woidPhasePortrait_speed_'...
        num2str(speed,'%1.0e') '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end
tilefigs()