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
revRatesClusterEdge = [0, 0.1, 0.2, 0.4, 0.8, 1.6];
speeds = [0.33];
slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125]);
% slowspeeds = fliplr([0.33, 0.2, 0.1, 0.05]);
attractionStrength = [0];
for speed = speeds
    phasePortraitFig = figure;
    plotCtr = 1;
    phasePortraitFig.Name = ['speed = ' num2str(speed)];
    for slowspeed = slowspeeds
        for revRateClusterEdge = revRatesClusterEdge
            filename = ['../results/woids/woids_N_' num2str(Nval) '_L_' num2str(Lval) ...
                '_noVolExcl'... %'_noUndulations'...
                '_v0_' num2str(speed,'%1.0e') ...
                '_vs_' num2str(slowspeed,'%1.0e') '_gradualSlowDown' ...
                '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e') '.mat'];
            if exist(filename,'file')
                load(filename)
                numFrames = size(xyarray,4);
                frames2plot = round(0.1*numFrames):numFrames; % discard first 10% as transient
                x = xyarray(:,:,1,frames2plot);
                y = xyarray(:,:,2,frames2plot);
                subplot(length(slowspeeds),length(revRatesClusterEdge),plotCtr)
                histogram2(x(:),y(:),'Normalization','Probability',...
                    'DisplayStyle','tile','EdgeColor','none','ShowEmptyBins','on')
                title(['r=' num2str(revRateClusterEdge) ', v_s =' num2str(slowspeed,'%1.0e')],...
                    'FontWeight','normal')
                ax = gca;
                ax.Position = ax.Position.*[1 1 1.2 1.2] - [0.0 0.0 0 0]; % stretch panel
                ax.DataAspectRatio = [1 1 1];
                ax.XTick = [];
                ax.YTick = [];
                ax.Box = 'on';
            end
            plotCtr = plotCtr + 1;
        end
    end
    %% export figure
    phasePortraitFig.PaperUnits = 'centimeters';
    filename = ['../figures/woids/woidPhasePortraitOccupancy_N_' num2str(Nval) '_L_' num2str(Lval)...
                '_noVolExcl'... %'_noUndulations'...
    '_speed_' num2str(speed,'%1.0e') '_slowing' '_gradual'...
        '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end