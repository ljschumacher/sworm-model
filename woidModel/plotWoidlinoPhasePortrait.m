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

M = 2;
radius = 0.35;
plotColor = [0.5, 0.5, 0.5];

revRatesClusterEdge = [0, 0.2, 0.4, 0.6, 0.8];
speeds = [0.15, 0.3];
attractionStrengths = [1e-3, 1e-4, 1e-5, 1e-6,0];
for speed = speeds
    phasePortraitFig = figure;
    plotCtr = 1;
    phasePortraitFig.Name = ['speed = ' num2str(speed)];
    for attractionStrength = attractionStrengths
        for revRateClusterEdge = revRatesClusterEdge
            filename = ['results/woidlinos/wlM' num2str(M) '_v0_' num2str(speed,'%1.0e') ...
                '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e') '_noContactForces.mat'];
            if exist(filename,'file')
                load(filename)
                positions2plot = xyarray(:,:,:,end);
                subplot(length(attractionStrengths),length(revRatesClusterEdge),plotCtr)
                ax = plotWoidTrajectoriesSingleFrame(positions2plot,L,radius,plotColor);
                title(['r=' num2str(revRateClusterEdge) ', \epsilon =' num2str(attractionStrength,'%1.0e')],...
                    'FontWeight','normal')
                ax.Position = ax.Position.*[1 1 1.2 1.2] - [0.0 0.0 0 0]; % stretch panel
                ax.DataAspectRatio = [1 1 1];
            end
            plotCtr = plotCtr + 1;
        end
    end
    %% export figure
    phasePortraitFig.PaperUnits = 'centimeters';
    filename = ['figures/woidlinos/woidlinoPhasePortrait_noContactForces_speed_'...
        num2str(speed,'%1.0e') '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end
tilefigs()