% plot woidlet phase portrait
% shows the end-points of of woidlet simulations for a 2D parameter sweep
close all
clear

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',14,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',9,...
    'LineWidth',2);

radius = 0.5;
plotColor = [0.5, 0.5, 0.5];

revRates = [0, 0.1, 1];
speeds = [0.1, 0.5, 1];
attractionStrengths = [1e-5, 1e-4, 5e-4];
for revRate = revRates
    phasePortraitFig = figure;
    plotCtr = 1;
    phasePortraitFig.Name = ['reversal rate = ' num2str(revRate)];
    for speed = speeds
        for attractionStrength = attractionStrengths
            filename = ['results/woidlets/wl_v0_' num2str(speed,'%1.0e') ...
                '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                '_revRate_' num2str(revRate,'%1.0e') '_noContactForces.mat'];
            if exist(filename,'file')
                load(filename)
                positions2plot = xyarray(:,:,:,end);
                subplot(length(speeds),length(attractionStrengths),plotCtr)
                ax = plotWoidTrajectoriesSingleFrame(positions2plot,L,radius,plotColor);
                title(['v=' num2str(speed) ', \epsilon =' num2str(attractionStrength,'%1.0e')],...
                    'FontWeight','normal')
                ax.Position = ax.Position.*[1 1 1.23 1.23] - [0.05 0.05 0 0]; % stretch panel
            end
            plotCtr = plotCtr + 1;
        end
    end
    %% export figure
    phasePortraitFig.PaperUnits = 'centimeters';
    filename = ['figures/woidlets/woidletPhasePortrait_noContactForces_revRate_'...
        num2str(revRate,'%1.0e') '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end
tilefigs()