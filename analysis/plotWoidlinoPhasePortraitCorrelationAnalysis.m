function [] = plotWoidlinoPhasePortraitCorrelationAnalysis(N,M,L,attractionStrength)
% plot woidlet phase portrait
% shows plots of correlation analysis

% issues/to-do:

close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',17,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',6,...
    'LineWidth',1,...
    'Renderer','opengl');

numRepeats = 1;
% revRatesClusterEdge = [0, 0.1, 0.2, 0.4, 0.8, 1.6];
drdN_rev_values = 0:0.2:1;
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125]);
% slowspeeds = fliplr([0.33, 0.05, 0.025, 0.0125]);
slowspeeds = [0.018];
trackedNodes = 1:max(round(M*0.16),1);
distBinwidth = 0.1; % in units of mm, sensibly to be chosen similar worm width or radius
maxDist = 2;
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = fliplr([0:0.2:1]);
% angleNoise = 1;

secondVariables = dkdN_dwell_values;

for speed = speeds
    poscorrFig = figure;
%     speedFig = figure;
%     dircorrFig = figure;
%     velcorrFig = figure;
    plotCtr = 1;
    poscorrFig.Name = ['v_0 = ' num2str(speed)];
    for slowspeed = slowspeeds
        for dkdN_dwell = dkdN_dwell_values
            for drdN_rev = drdN_rev_values
                gr = cell(numRepeats,1);
                for repCtr = 1:numRepeats
                    filename = ['../results/woidlinos/mapping/wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L) ...
                        ...'_noVolExcl' ...'_angleNoise_' num2str(angleNoise)...
                        '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                        '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
                        '_dkdN_' num2str(dkdN_dwell) ...
                        ...'_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                        '_revdensity_drdN_' num2str(drdN_rev) ...
                        '_run' num2str(repCtr) '.mat'];
                    if exist(filename,'file')
                        thisFile = load(filename);
                        maxNumFrames = size(thisFile.xyarray,4);
                        burnIn = round(500./thisFile.T*maxNumFrames);
                        if isfield(thisFile.param,'saveEvery')
                            saveEvery = thisFile.param.saveEvery;
                        else
                            saveEvery = thisFile.saveevery;
                        end
                        numFrames =  min(round((maxNumFrames - burnIn)*thisFile.param.dT*saveEvery/3),maxNumFrames - burnIn);
                        framesAnalyzed = burnIn + randperm(maxNumFrames - burnIn,numFrames); % randomly sample frames without replacement
                        %                 framesAnalyzed = round(linspace(burnIn,maxNumFrames,numFrames));
                        %                             framesAnalyzed = burnIn+1:maxNumFrames;
                        %% calculate stats
                        [~,~, ~,~, ~,~,...
                            gr{repCtr},distBins,~,~] = ...
                            correlationanalysisSimulations(thisFile,trackedNodes,distBinwidth,framesAnalyzed,maxDist);
%                         % plot lines for this file
%                         % speed v distance
%                         set(0,'CurrentFigure',speedFig)
%                         subplot(length(secondVariables),length(drdN_rev_values),plotCtr)
%                         boundedline(nearestDistBins,s_med,[s_med - s_ci(:,1), s_ci(:,2) - s_med],'alpha')
%                         if repCtr==1, hold on, end
%                         % directional and velocity cross-correlation
%                         set(0,'CurrentFigure',dircorrFig)
%                         subplot(length(secondVariables),length(drdN_rev_values),plotCtr)
%                         boundedline(pairDistBins,corr_o_med,[corr_o_med - corr_o_ci(:,1),...
%                             corr_o_ci(:,2) - corr_o_med],'alpha')
%                         if repCtr==1, hold on, end
%                         set(0,'CurrentFigure',velcorrFig)
%                         subplot(length(secondVariables),length(drdN_rev_values),plotCtr)
%                         boundedline(pairDistBins,corr_v_med,[corr_v_med - corr_v_ci(:,1),...
%                             corr_v_ci(:,2) - corr_v_med],'alpha')
                        if repCtr==1, hold on, end
                    end
                end
                %% combine data from multiple runs
                gr= horzcat(gr{:});
                %% plot data
                % radial distribution / pair correlation
                set(0,'CurrentFigure',poscorrFig)
                subplot(length(secondVariables),length(drdN_rev_values),plotCtr)
                boundedline(distBins(2:end)-distBinwidth/2,mean(gr,2),...
                    [std(gr,0,2) std(gr,0,2)])%./sqrt(size(gr,2)))
                ax = formatAxes(drdN_rev,dkdN_dwell);
                ax.YTick = 0:2:12;
                ax.YLim = [0 12];
%                 % format other plots
%                 % speed v distance
%                 set(0,'CurrentFigure',speedFig)
%                 subplot(length(secondVariables),length(drdN_rev_values),plotCtr)
%                 ax = formatAxes(drdN_rev,dkdN_dwell);
%                 ax.YLim = [0 0.5];
%                 ax.YTick = 0:0.1:0.5;
%                 ax.XDir = 'reverse';
%                 % directional and velocity cross-correlation
%                 set(0,'CurrentFigure',dircorrFig)
%                 subplot(length(secondVariables),length(drdN_rev_values),plotCtr)
%                 ax = formatAxes(drdN_rev,dkdN_dwell);
%                 ax.YLim = [-1 1];
%                 ax.YTick = [-1 0 1];
%                 set(0,'CurrentFigure',velcorrFig)
%                 subplot(length(secondVariables),length(drdN_rev_values),plotCtr)
%                 ax = formatAxes(drdN_rev,dkdN_dwell);
%                 ax.YLim = [-1 1];
%                 ax.YTick = [-1 0 1];
                % advance to next subplot
                plotCtr = plotCtr + 1;
            end
        end
    end
    %% export figures
    % radial distribution / pair correlation
    poscorrFig.PaperUnits = 'centimeters';
    fignameprefix = ['figures/woidlinoPhasePortrait'];
    fignamesuffix = ['_M' num2str(M)...
        'N_' num2str(thisFile.N) '_L_' num2str(thisFile.L(1)) ...
        '_noVolExcl' ...'_angleNoise_' num2str(angleNoise)...
        '_speed_' num2str(speed,'%1.0e') ...
        '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
        '_epsLJ_' num2str(attractionStrength,'%1.0e')...
        '.eps'];
    filename = [fignameprefix 'Radialdistribution' fignamesuffix];
    exportfig(poscorrFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
%     % speed v distance
%     speedFig.PaperUnits = 'centimeters';
%     filename = [fignameprefix 'SpeedvDistance' fignamesuffix];
%     exportfig(speedFig,filename, exportOptions)
%     system(['epstopdf ' filename]);
%     system(['rm ' filename]);
%     % directional and velocity cross-correlation
%     dircorrFig.PaperUnits = 'centimeters';
%     filename = [fignameprefix 'Dirxcorr' fignamesuffix];
%     exportfig(dircorrFig,filename, exportOptions)
%     system(['epstopdf ' filename]);
%     system(['rm ' filename]);
%     velcorrFig.PaperUnits = 'centimeters';
%     filename = [fignameprefix 'Velxcorr' fignamesuffix];
%     exportfig(velcorrFig,filename, exportOptions)
%     system(['epstopdf ' filename]);
%     system(['rm ' filename]);
end
end

function ax = formatAxes(revRateClusterEdge,var2)
title(['r=' num2str(revRateClusterEdge) ', dk/dN =' num2str(var2)],...
    'FontWeight','normal')
ax = gca;
ax.Position = ax.Position.*[1 1 1.2 1.2] - [0.0 0.0 0 0]; % stretch panel
ax.XLim = [0 2];
ax.XTick = 0:2;
ax.XTickLabel = [];
ax.YTickLabel = [];
ax.Box = 'on';
grid on
end