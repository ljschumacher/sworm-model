% plot correlation analysis for simulations


clear
close all

% set pdf export options
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1);

% set plotting options
distBinWidth = 0.1; % in units of mm, sensibly to be chosen similar worm width or radius
maxDist = 2;
plotbins = (0:distBinWidth:(maxDist-distBinWidth)) + distBinWidth/2;
nBins = numel(plotbins);
ylabels = {'directional cross-corr.','velocity cross-corr.','vel. to nbr dir. corr.'};
figurepath = 'figures/';
figurePrefices = {'dircrosscorr/dircrosscorr_','velcrosscorr/velcrosscorr_','velnbrcorr/velnbrcorr_'};

% set simulation file
filepath = '~/Dropbox/projects/collectiveBehaviour/sworm-model/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r2/postiPredictiveCheck/';
filenames = rdir([filepath '*run1.mat']);
nFiles = numel(filenames);
numReps = 5;

% set analysis parameters
trackedNodes = 1:3;

% loop over file to plot
for fileCtr = 1:nFiles
    fileName = strrep(filenames(fileCtr).name,filepath,'');
    corr_o_mean = NaN(numReps,nBins);
    corr_o_ci = NaN(numReps,nBins,2);
    corr_v_mean = NaN(numReps,nBins);
    corr_vn_mean = NaN(numReps,nBins);
    corr_v_ci = NaN(numReps,nBins,2);
    corr_vn_ci = NaN(numReps,nBins,2);
    if exist([filepath fileName],'file')
        for repCtr = 1:numReps
            thisFileName = strrep(fileName,'run1',['run' num2str(repCtr)]);
            thisFile = load([filepath thisFileName]);
            %% analyse simulation data
            % pick frames to analyse, e.g. second half of simulation
            maxNumFrames = size(thisFile.xyarray,4);
            burnIn = round(0.5*maxNumFrames);
            framesAnalyzed = burnIn+1:maxNumFrames;
            % calculate stats
            [corr_o_mean(repCtr,:),corr_o_ci(repCtr,:,:),...
                corr_v_mean(repCtr,:),corr_v_ci(repCtr,:,:),...
                corr_vn_mean(repCtr,:),corr_vn_ci(repCtr,:,:),~,~,nbrDistBins,pairDistBins] = ...
                correlationanalysisSimulations(thisFile,trackedNodes,distBinWidth,framesAnalyzed,maxDist);
        end
        %% plot analysis results
        dircorrFig = figure; hold on
        velcorrFig = figure; hold on
        velncorrFig = figure; hold on
        % directional and velocity cross-correlation
        boundedline(plotbins,mean(corr_o_mean),...
            [mean(corr_o_mean - corr_o_ci(:,:,1)); mean(corr_o_ci(:,:,2) - corr_o_mean)]',...
            dircorrFig.Children)
        boundedline(plotbins,mean(corr_v_mean),...
            [mean(corr_v_mean - corr_v_ci(:,:,1)); mean(corr_v_ci(:,:,2) - corr_v_mean)]',...
            velcorrFig.Children)
        boundedline(plotbins,mean(corr_vn_mean),...
            [mean(corr_vn_mean - corr_vn_ci(:,:,1)); mean(corr_vn_ci(:,:,2) - corr_vn_mean)]',...
            velncorrFig.Children)
        %% format and export figures
        figHandles = [dircorrFig, velcorrFig, velncorrFig];
        for  figCtr = 1:length(figHandles)% common formating for all figures
            figHandle = figHandles(figCtr);
            figHandle.Children.XLim = [0 maxDist];
            set(figHandle,'PaperUnits','centimeters')
            plot(figHandle.Children,plotbins,zeros(size(plotbins)),'k--')
            xlabel(figHandle.Children,'distance r (mm)')
            ylabel(figHandle.Children,ylabels{figCtr})
            % export figure
            figurename = [figurepath figurePrefices{figCtr} ...
                strrep(strrep(fileName,'run1',''),'.mat','')];
            exportfig(figHandle,[figurename '.eps'],exportOptions)
            system(['epstopdf ' figurename '.eps']);
            system(['rm ' figurename '.eps']);
        end
    end
end