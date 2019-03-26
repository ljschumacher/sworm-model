% plot correlation analysis for simulations


clear
close all

% set pdf export options
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

% set plotting options
nBins = 10;
figurepath = 'figures/speedratio/';

% set simulation file
filepath = '~/Dropbox/projects/collectiveBehaviour/sworm-model/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r2/postiPredictiveCheck/';
filenames = rdir([filepath '*M18*v0_0.14*run1.mat']);
nFiles = numel(filenames);
numReps = 10;

% set analysis parameters
trackedNodes = 1:3;

% loop over file to plot
for fileCtr = 1:nFiles
    fileName = strrep(filenames(fileCtr).name,filepath,'');
    speedRatios = NaN(numReps,nBins,2);
    if exist([filepath fileName],'file')
        for repCtr = 1:numReps
            thisFileName = strrep(fileName,'run1',['run' num2str(repCtr)]);
            thisFile = load([filepath thisFileName]);
            %% analyse simulation data
            % pick frames to analyse, e.g. second half of simulation
            maxNumFrames = size(thisFile.xyarray,4);
            burnIn = round(0.5*maxNumFrames);
            framesAnalyzed = (burnIn+1):maxNumFrames;
            % calculate stats
            [speedRatios(repCtr,:,:),densitybinedges] = ...
                speedratioAnalysisSimulations(thisFile,trackedNodes,framesAnalyzed,nBins);
        end
        %% plot analysis results
        figure; hold on
        errorbar(densitybinedges(1:nBins),nanmean(speedRatios(:,:,2))./nanmean(speedRatios(:,:,1)),...
            sqrt((nanstd(speedRatios(:,:,2))./nanmean(speedRatios(:,:,2))).^2 ...
            + (nanstd(speedRatios(:,:,1))./nanmean(speedRatios(:,:,1))).^2) ...
        .*nanmean(speedRatios(:,:,2))./nanmean(speedRatios(:,:,1))./sqrt(numReps))
        
        %% format and export figures
        figHandle = gcf;
        box(figHandle.Children,'on')
        set(figHandle,'PaperUnits','centimeters')
        ylabel('P_f/P_s')
        xlabel('local density, \rho_6 (worms/mm^2)')
%         ylim([0 ceil(10*max(nanmean(speedRatios(:,:,2))./nanmean(speedRatios(:,:,1))))/10])
        % export figure
        figurename = [figurepath 'speedRatio_' ...
            strrep(strrep(fileName,'run1',''),'.mat','')];
        exportfig(figHandle,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
    end
end