% plot correlation analysis

% issues/to-do:

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1,...
    'Renderer','opengl');

distBinwidth = 0.035; % in units of mm, sensibly to be chosen similar worm width or radius
maxDist = 2;
% simulations = {'DA609_noflux','N2_noflux','DA609_noflux_slowingNodesAll',...
%     'DA609_noflux_lennardjones1e-04'};
resultsfiles = rdir('../results/woids/*N_40_L_7.5_v0*.mat');
for ii = 1:length(resultsfiles)
    simulations{ii} = strrep(strrep(resultsfiles(ii).name,'../results/woids/woids_',''),'.mat','');
end
nSims = length(simulations);
plotColors = lines(nSims);
trackedNodesNames = {'head'};%,'body'};
trackedNodesDict = containers.Map({'head','body'},{1:8; 1:49});% which nodes to calculate the tracking stats from, to compare eg with pharynx labeled expmntal data

for trackedNodesName = trackedNodesNames
    trackedNodes = trackedNodesDict(trackedNodesName{1});
    for simCtr = 1:nSims % can be parfor
        filenames = {[simulations{simCtr}]};
        speedFig = figure; hold on
        dircorrFig = figure; hold on
        velcorrFig = figure; hold on
        poscorrFig = figure; hold on
        %% load data
        numFiles = length(filenames);
        for fileCtr = 1:numFiles
% %             filename = ['../results/' filenames{fileCtr} '.mat'];
            filename = ['../results/woids/woids_' filenames{fileCtr} '.mat'];
            thisFile = load(filename);
            maxNumFrames = size(thisFile.xyarray,4);
            burnIn = round(0.1*maxNumFrames);
            numFrames =  min(round((maxNumFrames - burnIn)*thisFile.param.dT*thisFile.saveevery),maxNumFrames - burnIn); %maxNumFrames-burnIn;
            framesAnalyzed = burnIn + randperm(maxNumFrames - burnIn,numFrames); % randomly sample frames without replacement
%             framesAnalyzed = round(linspace(burnIn,maxNumFrames,numFrames));
%             framesAnalyzed = burnIn+1:maxNumFrames;
            %% calculate stats
            [s_med,s_ci, corr_o_med,corr_o_ci, corr_v_med,corr_v_ci, gr,distBins,nearestDistBins,pairDistBins] = ...
    correlationanalysisSimulations(thisFile,trackedNodes,distBinwidth,framesAnalyzed,maxDist);
            %% plot data
            % speed v distance
            boundedline(nearestDistBins,s_med,[s_med - s_ci(:,1), s_ci(:,2) - s_med],...
                'alpha',speedFig.Children,'cmap',plotColors(simCtr,:))
            % directional and velocity cross-correlation
            boundedline(pairDistBins,corr_o_med,[corr_o_med - corr_o_ci(:,1), corr_o_ci(:,2) - corr_o_med],...
                'alpha',dircorrFig.Children,'cmap',plotColors(simCtr,:))
            boundedline(pairDistBins,corr_v_med,[corr_v_med - corr_v_ci(:,1), corr_v_ci(:,2) - corr_v_med],...
                'alpha',velcorrFig.Children,'cmap',plotColors(simCtr,:))
            % radial distribution / pair correlation
            boundedline(distBins(2:end)-distBinwidth/2,mean(gr,2),...
                [nanstd(gr,0,2) nanstd(gr,0,2)]./sqrt(numFrames),...
                'alpha',poscorrFig.Children,'cmap',plotColors(simCtr,:))
        end
        %% format and export figures
        for figHandle = [speedFig, dircorrFig, velcorrFig, poscorrFig] % common formating for both figures
            title(figHandle.Children,{simulations{simCtr} ; [trackedNodesName{1} ' tracked']},...
                'FontWeight','normal','Interpreter','none');
            set(figHandle,'PaperUnits','centimeters')
        end
        % speed v distance
        speedFig.Children.XLim = [0 2];
        ylabel(speedFig.Children,'speed (mm/s)')
        xlabel(speedFig.Children,'distance to nearest neighbour (mm)')
        speedFig.Children.Box = 'on';
        speedFig.Children.XDir = 'reverse';
        figurename = ['figures/speedvdistance/speedvsneighbourdistance_' simulations{simCtr} '_' trackedNodesName{1}];
        exportfig(speedFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        % directional and velocity cross-correlation
%         dircorrFig.Children.YLim = [-1 1];
        dircorrFig.Children.XLim = [0 2];
        ylabel(dircorrFig.Children,'directional correlation')
        xlabel(dircorrFig.Children,'distance r (mm)')
        figurename = ['figures/dircrosscorr/dircrosscorr' simulations{simCtr} '_' trackedNodesName{1}];
        exportfig(dircorrFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        %
%         velcorrFig.Children.YLim = [-1 1];
        velcorrFig.Children.XLim = [0 2];
        ylabel(velcorrFig.Children,'velocity correlation')
        xlabel(velcorrFig.Children,'distance r (mm)')
        figurename = ['figures/velcrosscorr/velcrosscorr' simulations{simCtr} '_' trackedNodesName{1}];
        exportfig(velcorrFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        % radial distribution / pair correlation
        poscorrFig.Children.XLim = [0 2];
        ylabel(poscorrFig.Children,'positional correlation g(r)')
        xlabel(poscorrFig.Children,'distance r (mm)')
        figurename = ['figures/radialdistribution/radialdistributionfunction_' simulations{simCtr} '_' trackedNodesName{1}];
        exportfig(poscorrFig,[figurename '.eps'],exportOptions)
        system(['epstopdf ' figurename '.eps']);
        system(['rm ' figurename '.eps']);
        %%
        close all
    end
end