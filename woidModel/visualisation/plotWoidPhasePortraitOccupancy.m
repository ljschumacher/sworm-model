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
    'FontSize',6,...
    'LineWidth',1,...
    'Renderer','opengl');

centering = true;

radius = 0.035;
plotColor = [0.25, 0.25, 0.25];
Nval = 40;
Lval = 7.5;
revRatesClusterEdge = [0, 0.4, 0.8, 1.6, 3.2];
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125, 0.005]);
% slowspeeds = fliplr([0.33, 0.025, 0.0125, 0.005, 0.001]);
slowspeeds = [0.018];
attractionStrength = [0];
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = fliplr([0 1./[8 4 2 1]]);
secondVariables = dkdN_dwell_values;
for speed = speeds
    phasePortraitFig = figure;
    plotCtr = 1;
    phasePortraitFig.Name = ['speed = ' num2str(speed)];
    for slowspeed = slowspeeds
        for dkdN_dwell = dkdN_dwell_values
        for revRateClusterEdge = revRatesClusterEdge
                filename = ['../results/woids/woids_N_' num2str(Nval) '_L_' num2str(Lval) ...
                    ...'_noUndulations'...'_noVolExcl' ...'_angleNoise'
                    '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e')...
                    '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
                    '_dkdN_' num2str(dkdN_dwell) ...num2str(num_nbr_max_per_nodes)...
                    '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                    '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e')...
                    '_run1.mat'];
            if exist(filename,'file')
                load(filename)
                numFrames = size(xyarray,4);
                frames2plot = numFrames-round(250./T*numFrames):numFrames; % plots last 250 seconds
                if centering&&numel(L)==2
                    xyarray(:,:,1,frames2plot) = xyarray(:,:,1,frames2plot) + L(1)/2 - mean(mean(mean(xyarray(:,:,1,frames2plot),4),2),1);
                    xyarray(:,:,2,frames2plot) = xyarray(:,:,2,frames2plot) + L(2)/2 - mean(mean(mean(xyarray(:,:,2,frames2plot),4),2),1);
                    for frameCtr = frames2plot
                        [ xyarray(:,:,:,frameCtr), ~ ] = checkWoidBoundaryConditions(xyarray(:,:,:,frameCtr), [], 'periodic', L);
                    end
                end
                x = xyarray(:,:,1,frames2plot);
                y = xyarray(:,:,2,frames2plot);
                subplot(length(secondVariables),length(revRatesClusterEdge),plotCtr)
                histogram2(x(:),y(:),'Normalization','Probability',...
                    'DisplayStyle','tile','EdgeColor','none','ShowEmptyBins','on')
                title(['r=' num2str(revRateClusterEdge) ', dk/dN =' num2str(dkdN_dwell,'%1.0e')],...
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
    end
    %% export figure
    phasePortraitFig.PaperUnits = 'centimeters';
    filename = ['../figures/woids/woidPhasePortraitOccupancy_N_' num2str(Nval) '_L_' num2str(Lval)...
        ...'_noUndulations'...'_noVolExcl' ...'_angleNoise'
        '_speed_' num2str(speed,'%1.0e') '_slowing' '_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
        '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end