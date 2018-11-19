% plot woidlino phase portrait
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

M = 18;
L = [3.6 3.6];
N = 40;
plotColor = [0.25, 0.25, 0.25];

revRatesClusterEdge = 0:5;

speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.05, 0.025, 0.0125]);
slowspeeds = [0.018];
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = fliplr(0:0.2:1);
dkdN_undwell_values = fliplr(0:0.2:2);
% angleNoise = 1;
k_theta = 2;
% f_hapt = 0.5;

secondVariables = dkdN_dwell_values;
ndwellVals = numel(dkdN_dwell_values);
nundwellVals = numel(dkdN_undwell_values);
aspectRatio = nundwellVals/(ndwellVals);

for speed = speeds
    for revRateClusterEdge = revRatesClusterEdge
        
        phasePortraitFig = figure;
        plotCtr = 1;
        phasePortraitFig.Name = ['speed = ' num2str(speed)];
        for slowspeed = slowspeeds
            for dkdN_dwell = dkdN_dwell_values
                for dkdN_undwell = dkdN_undwell_values
                    filename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
                        ...'_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
                        '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                        '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
                        '_dkdN_' num2str(dkdN_dwell) '_' num2str(dkdN_undwell)... %% also check the case of indem dwell/roam
                        '_revRateClusterEdge_' num2str(revRateClusterEdge,2) ...
                        ...'_haptotaxis_' num2str(f_hapt) ...
                        '_clusteredStart' ...
                        '_run1.mat'];
                    filepath = '../results/woidlinos/';
                    if exist([filepath filename],'file')
                        load([filepath filename])
                        time2plot = round(size(xyarray,4)*(0.9 + 0.1*rand()));
                        positions2plot = xyarray(:,:,:,time2plot);
                        subplot(length(secondVariables),length(dkdN_undwell_values),plotCtr)
                        radius = param.rc;
                        Rgyr = sqrt(sum(var(positions2plot(:,1,:))));
                        thisColor = [1 0 0]*log10(Rgyr)/2 + lines(1)*(1 - log10(Rgyr)/2);
                        ax = plotWoidTrajectoriesSingleFrame(positions2plot,[],radius,thisColor,true);
                        ax.Position = ax.Position.*[1 1 1.25 1.25] - [0.0 0.0 0 0]; % stretch panel
                        %                         ax.DataAspectRatio = [1 1 1];
                        ax.Box = 'on';
                    else
                        warning([filename ' does not exist'])
                    end
                    plotCtr = plotCtr + 1;
                end
            end
        end
        % make overall axes
        ax = axes('Color','none');
        ax.XTick = linspace(1/nundwellVals/2,1-1/nundwellVals/2,nundwellVals);
        ax.XTickLabel = num2str(fliplr(dkdN_undwell_values)');
        ax.YTick = linspace(1/ndwellVals/2,1-1/ndwellVals/2,ndwellVals);
        ax.YTickLabel = num2str(fliplr(dkdN_dwell_values)');
        ax.TickDir = 'out';
        ax.Position = ax.Position.*[1 1 1 1];
        xlabel('dk_{roam}/d\rho')
        ylabel('dk_{dwell}/d\rho')
        %% export figure
        phasePortraitFig.Position(3) = phasePortraitFig.Position(4)*aspectRatio; % resize figure
        phasePortraitFig.PaperUnits = 'centimeters';
        filename = ['../figures/woidlinos/woidlinoPhaseDiagram_clusterStability_N_' num2str(N) ...
            '_M_' num2str(M) '_L_' num2str(L(1)) '_noVolExcl' ...
            ...'_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
            '_speed_' num2str(speed,'%1.0e') ...
            '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
            '_independent_revRate_' num2str(revRateClusterEdge) ...'_haptotaxis_' num2str(f_hapt) ...
            '.eps'];
        exportfig(phasePortraitFig,filename, exportOptions)
        system(['epstopdf ' filename]);
        system(['rm ' filename]);
    end
end
