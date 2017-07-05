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
    'Renderer','painters');

M = 2;
L = [7.5 7.5];
N = 40;
plotColor = [0.5, 0.5, 0.5];

revRatesClusterEdge = [0, 0.1, 0.2, 0.4, 0.8];
speeds = [0.33];
slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025]);
attractionStrengths = [0];
for speed = speeds
    phasePortraitFig = figure;
    plotCtr = 1;
    phasePortraitFig.Name = ['speed = ' num2str(speed)];
    for slowspeed = slowspeeds
        for attractionStrength = attractionStrengths
            for revRateClusterEdge = revRatesClusterEdge
                filename = ['../results/woidlinos/wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
                    '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                    '_gradualSlowDown' ...
                    '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                    '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e') '.mat'];
                if exist(filename,'file')
                    load(filename)
                    time2plot = round(size(xyarray,4)*(0.9 + 0.1*rand()));
                    positions2plot = xyarray(:,:,:,time2plot);
                    subplot(length(slowspeeds),length(revRatesClusterEdge),plotCtr)
                    radius = param.rc;
                    ax = plotWoidTrajectoriesSingleFrame(positions2plot,L,radius,plotColor);
                    title(['r=' num2str(revRateClusterEdge) ', v_s =' num2str(slowspeed,'%1.0e')],...
                        'FontWeight','normal')
                    ax.Position = ax.Position.*[1 1 1.2 1.2] - [0.0 0.0 0 0]; % stretch panel
                    ax.DataAspectRatio = [1 1 1];
                    ax.Box = 'on';
                else
                    warning([filename ' does not exist'])
                end
                plotCtr = plotCtr + 1;
            end
        end
    end
    %% export figure
    phasePortraitFig.PaperUnits = 'centimeters';
    filename = ['../figures/woidlinos/woidlinoPhasePortrait_N_' num2str(N) ...
        '_M_' num2str(M) '_L_' num2str(L(1))...
        '_speed_' num2str(speed,'%1.0e') '_slowing' '_gradual' ...
        '.eps'];
    exportfig(phasePortraitFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end
