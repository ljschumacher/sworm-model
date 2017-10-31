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
revRatesClusterEdge = [0, 0.2, 0.4, 0.8, 1.6, 3.2];
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125]);
% slowspeeds = fliplr([0.33, 0.025, 0.0125, 0.005, 0.001]);
slowspeeds = [0.018];
attractionStrength = [0];
slowingMode = 'stochastic';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = fliplr([0 1./[8 4 2 1 0.5]]);
% num_nbr_max_per_nodes = 4;
for speed = speeds
    phasePortraitFig = figure;
    plotCtr = 1;
    phasePortraitFig.Name = ['speed = ' num2str(speed)];
    for slowspeed = slowspeeds
        for dkdN_dwell = dkdN_dwell_values
            for revRateClusterEdge = revRatesClusterEdge
                filename = ['../results/woids/woids_N_' num2str(Nval) '_L_' num2str(Lval) ...
                    ...'_noVolExcl' ...'_angleNoise' ...%'_noUndulations'
                    '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e')...
                    '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
                    '_dkdN_' num2str(dkdN_dwell) ...num2str(num_nbr_max_per_nodes)...
                    '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                    '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e')...
                    '.mat'];
                if exist(filename,'file')
                    load(filename)
                    time2plot = round(size(xyarray,4)*(0.9 + 0.1*rand()));
                    positions2plot = xyarray(:,:,:,time2plot);
                    subplot(length(dkdN_dwell_values),length(revRatesClusterEdge),plotCtr)
                    ax = plotWoidTrajectoriesSingleFrame(positions2plot,L,radius,plotColor);
                    title(['r=' num2str(revRateClusterEdge) ', dk/dN =' num2str(dkdN_dwell,'%1.0e')],...
                        'FontWeight','normal')
                    ax.Position = ax.Position.*[1 1 1.2 1.2] - [0.0 0.0 0 0]; % stretch panel
                    ax.DataAspectRatio = [1 1 1];
                    ax.Box = 'on';
                end
                plotCtr = plotCtr + 1;
            end
        end
    end
    %% export figure
    phasePortraitFig.PaperUnits = 'centimeters';
    filename = ['../figures/woids/woidPhasePortrait_N_' num2str(Nval) '_L_' num2str(Lval) ...
        ...'_noVolExcl' ...'_angleNoise' ...%'_noUndulations'
        '_speed_' num2str(speed,'%1.0e') '_slowing' '_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...num2str(num_nbr_max_per_nodes) ...
        '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end