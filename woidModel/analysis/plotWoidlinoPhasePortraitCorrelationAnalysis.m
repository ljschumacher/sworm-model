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
revRatesClusterEdge = [0, 0.1, 0.2, 0.4, 0.8, 1.6];
speeds = [0.33];
slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125]);
trackedNodes = 1;
distBinwidth = 0.0525; % in units of mm, sensibly to be chosen similar worm width or radius
maxDist = 2;

for speed = speeds
    poscorrFig = figure;
    plotCtr = 1;
    poscorrFig.Name = ['v_0 = ' num2str(speed)];
    for slowspeed = slowspeeds
        for revRateClusterEdge = revRatesClusterEdge
            gr = cell(numRepeats,1);
            for repCtr = 1:numRepeats
                filename = ['../results/woidlinos/wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L) ...
                    '_noVolExcl' ...    
                    '_v0_' num2str(speed,'%1.0e') ...
                    '_vs_' num2str(slowspeed,'%1.0e') '_gradualSlowDown' ...
                    '_epsLJ_' num2str(attractionStrength,'%1.0e') ...
                    '_revRateClusterEdge_' num2str(revRateClusterEdge,'%1.0e') ...
                    '_run' num2str(repCtr) '.mat'];
                if exist(filename,'file')
                    thisFile = load(filename);
                    maxNumFrames = size(thisFile.xyarray,4);
                    burnIn = round(0.5*maxNumFrames);
                    if isfield(thisFile.param,'saveEvery')
                        saveEvery = thisFile.param.saveEvery;
                    else
                        saveEvery = thisFile.saveevery;
                    end
                    numFrames =  min(round((maxNumFrames - burnIn)*thisFile.param.dT*saveEvery),maxNumFrames - burnIn); %maxNumFrames-burnIn;
                    framesAnalyzed = burnIn + randperm(maxNumFrames - burnIn,numFrames); % randomly sample frames without replacement
                    %                 framesAnalyzed = round(linspace(burnIn,maxNumFrames,numFrames));
                    %                             framesAnalyzed = burnIn+1:maxNumFrames;
                    %% calculate stats
                    [~,~, ~,~, ~,~, gr{repCtr},distBins,~,~] = ...
                        correlationanalysisSimulations(thisFile,trackedNodes,distBinwidth,framesAnalyzed,maxDist);
                end
            end
            %% combine data from multiple runs
            gr= horzcat(gr{:});
            %% plot data
            % radial distribution / pair correlation
            set(0,'CurrentFigure',poscorrFig)
            subplot(length(slowspeeds),length(revRatesClusterEdge),plotCtr)
            boundedline(distBins(2:end)-distBinwidth/2,mean(gr,2),...
                [std(gr,0,2) std(gr,0,2)]./sqrt(size(gr,2)))
            ax = formatAxes(revRateClusterEdge,slowspeed);
            ax.YTick = 0:2:12;
            ax.YLim = [0 12];
            plotCtr = plotCtr + 1;
        end
    end
    %% export figures
    % radial distribution / pair correlation
    poscorrFig.PaperUnits = 'centimeters';
    filename = ['figures/woidlinoPhasePortraitRadialdistribution_M' num2str(M)...
        'N_' num2str(thisFile.N) '_L_' num2str(thisFile.L(1)) ...
        '_noVolExcl' ...%'_noUndulations'...
        '_speed_'...
        num2str(speed,'%1.0e') '_slowing' '_gradual'...
        '_epsLJ_' num2str(attractionStrength,'%1.0e')...
        '.eps'];
    exportfig(poscorrFig,filename, exportOptions)
    system(['epstopdf ' filename]);
    system(['rm ' filename]);
end
end

function ax = formatAxes(revRateClusterEdge,slowspeed)
title(['r=' num2str(revRateClusterEdge) ', v\_s =' num2str(slowspeed,'%1.0e')],...
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