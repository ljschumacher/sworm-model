% calculate speed vs neighbour distance, directional correlation, and
% radial distribution functions

% issues/to-do:
% - seperate into individual functions for each statistic?
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

% define functions for grpstats
mad1 = @(x) mad(x,1); % median absolute deviation
distBinwidth = 0.055; % in units of mm
simulations = {'DA609_noflux','DA609_noflux_slowingNodesAll','N2_noflux'};
trackedNodesList = {1:5; 1:50}; % which nodes to calculate the tracking stats from, to compare eg with pharynx labeled expmntal data
trackedNodesNames = {'head','body'};
for nodesCtr = 1:length(trackedNodesList)
    trackedNodes = trackedNodesList{nodesCtr};
for simCtr = 1:length(simulations)
filenames = {[simulations{simCtr}]};
speedFig = figure; hold on
dircorrFig = figure; hold on
poscorrFig = figure; hold on
%% load data
numFiles = length(filenames);
for fileCtr = 1:numFiles
    filename = ['../results/' filenames{fileCtr} '.mat'];
    load(filename);
    maxNumFrames = size(xyphiarray,4);
    numFrames = 1000;%round(maxNumFrames/4)%*dT);
    burnIn = round(0.1*maxNumFrames);
    framesAnalyzed = burnIn + randperm(maxNumFrames - 1 - burnIn,numFrames); % randomly sample frames without replacement
    
    %% calculate stats
    x = squeeze(mean(xyphiarray(:,trackedNodes,1,:),2)); % centroid of worm head
    y = squeeze(mean(xyphiarray(:,trackedNodes,2,:),2)); % centroid of worm head
    u = (x(:,2:end) - x(:,1:end-1))/dT;
    v = (y(:,2:end) - y(:,1:end-1))/dT;
    speed = sqrt(u(:,framesAnalyzed).^2 + v(:,framesAnalyzed).^2);
    N = size(x,1);
    pairdist = NaN(N*(N-1)/2,numFrames);
    dxcorr = NaN(size(pairdist));
    distBins = 0:distBinwidth:L;
    gr = NaN(length(distBins) - 1,numFrames);
    mindist = NaN(N,numFrames);
    for frameCtr = 1:numFrames
        frame = framesAnalyzed(frameCtr);
        dxcorr(:,frameCtr) = vectorCrossCorrelation2D(u(:,frame),v(:,frame),true); % directional correlation
        pairdist(:,frameCtr) = pdist([x(:,frame) y(:,frame)]); % distance between all pairs, in micrometer
        gr(:,frameCtr) = histcounts(pairdist(:,frameCtr),distBins,'Normalization','probability'); % radial distribution function
        gr(:,frameCtr) = gr(:,frameCtr)'.*L^2./(2*distBins(2:end)*distBinwidth); % normalization
        D = squareform(pairdist(:,frameCtr)); % distance of every worm to every other
        mindist(:,frameCtr) = min(D + max(max(D))*eye(size(D)));
    end
    [s_med,s_mad] = grpstats(speed(:),quant(mindist(:),distBinwidth),{@median,mad1});
    [c_med,c_mad] = grpstats(dxcorr(:),quant(pairdist(:),distBinwidth),{@median,mad1});
    %% plot data
    bins = (0:numel(s_med)-1).*distBinwidth;
    boundedline(bins,s_med,[s_mad, s_mad],...
        'alpha','transparency',1/max(numFiles,2),speedFig.Children)
    bins = (0:numel(c_med)-1).*distBinwidth;
    boundedline(bins,c_med,[c_mad, c_mad],...
        'alpha','transparency',1/max(numFiles,2),dircorrFig.Children)
    boundedline(distBins(2:end)-distBinwidth/2,mean(gr,2),[mad(gr,0,2) mad(gr,0,2)],...
        'alpha','transparency',1/max(numFiles,2),poscorrFig.Children)
end
%% format and export figures
for figHandle = [speedFig, dircorrFig, poscorrFig] % common formating for both figures
    title(figHandle.Children,[simulations{simCtr} ' simulation, ' trackedNodesNames{nodesCtr} ' tracked'],...
        'FontWeight','normal','Interpreter','none');
    set(figHandle,'PaperUnits','centimeters')
end
%
speedFig.Children.YLim = [0 0.4];
speedFig.Children.XLim = [0 2];
ylabel(speedFig.Children,'speed (mm/s)')
xlabel(speedFig.Children,'distance to nearest neighbour (mm)')
figurename = [simulations{simCtr} '_' trackedNodesNames{nodesCtr} '_speedvsneighbourdistance'];
exportfig(speedFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);
%
dircorrFig.Children.YLim = [-1 1];
dircorrFig.Children.XLim = [0 2];
ylabel(dircorrFig.Children,'directional cross-correlation')
xlabel(dircorrFig.Children,'distance between pair (mm)')
figurename = [simulations{simCtr} '_' trackedNodesNames{nodesCtr} '_dircrosscorr'];
exportfig(dircorrFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);
%
poscorrFig.Children.YLim = [0 7];
poscorrFig.Children.XLim = [0 2];
ylabel(poscorrFig.Children,'radial distribution function g(r)')
xlabel(poscorrFig.Children,'distance r (mm)')
figurename = [simulations{simCtr} '_' trackedNodesNames{nodesCtr} '_radialdistributionfunction'];
exportfig(poscorrFig,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);
end
end