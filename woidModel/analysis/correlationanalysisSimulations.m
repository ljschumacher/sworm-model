function [s_med,s_ci, corr_o_med,corr_o_ci, corr_v_med,corr_v_ci, gr,distBins,nNbrDistBins,pairDistBins] = ...
    correlationanalysisSimulations(simfile,trackedNodes,distBinWidth,framesAnalyzed,maxDist)
% calculate directional & velocity correlation, and
% radial distribution functions

% issues/to-do:
% - for periodic bounadry conditions, if the tracked nodes overlap a
% boundary, it's not clear how to best find the centroid
% - Area in g(r) needs to be adjusted if arena is circular

% define functions for grpstats
mad1 = @(x) mad(x,1); % median absolute deviation
% alternatively could use boxplot-style confidence intervals on the mean,
% which are 1.57*iqr/sqrt(n) - unclear how justified this is
iqrci = @(x) 1.57*iqr(x)/sqrt(numel(x));
% or one could use a bootstrapped confidence interval
bootserr = @(x) bootci(1e2,{@nanmedian,x},'alpha',0.05,'Options',struct('UseParallel',false));

% convert result to double precision
simfile.xyarray = double(simfile.xyarray);

M = size(simfile.xyarray,2);
if nargin<2||isempty(trackedNodes)
    trackedNodes = 1:M;
end
if nargin<3||isempty(distBinWidth)
    distBinWidth = 0.05; % in units of mm
end
if nargin<4||isempty(framesAnalyzed)
    framesAnalyzed = 1:size(simfile.xyarray,4);
end
if nargin<5||isempty(maxDist)
    maxDist = min(simfile.L)/2;
end
numFrames = numel(framesAnalyzed);
if isfield(simfile.param,'saveEvery')
    saveEvery = simfile.param.saveEvery;
else
    saveEvery = simfile.saveevery;
end
%% calculate stats
if strcmp(simfile.param.bc,'periodic')
    x = squeeze(simfile.xyarray(:,round(mean(trackedNodes)),1,:)); % centroid of tracked obj
    y = squeeze(simfile.xyarray(:,round(mean(trackedNodes)),2,:)); % centroid of tracked obj
else
    x = squeeze(mean(simfile.xyarray(:,trackedNodes,1,:),2)); % centroid of tracked obj
    y = squeeze(mean(simfile.xyarray(:,trackedNodes,2,:),2)); % centroid of tracked obj
end
dxdt = diff(x,1,2);
dydt = diff(y,1,2);
if numel(trackedNodes)>1
    dxds = diff(simfile.xyarray(:,fliplr(trackedNodes),1,:),1,2);
    dyds = diff(simfile.xyarray(:,fliplr(trackedNodes),2,:),1,2);
    if strcmp(simfile.param.bc,'periodic')
        dxdt = correctForPeriodicBoundary(dxdt,simfile.L(1));
        dydt = correctForPeriodicBoundary(dydt,simfile.L(2));
        dxds = correctForPeriodicBoundary(dxds,simfile.L(1));
        dyds = correctForPeriodicBoundary(dyds,simfile.L(2));
    end
    ox = squeeze(mean(dxds,2));
    oy = squeeze(mean(dyds,2));
end
% central difference (should be equibalent to matlab's gradient function
dxdt = ([dxdt, dxdt(:,end)] + [dxdt(:,1), dxdt])./2;
dydt = ([dydt, dydt(:,end)] + [dydt(:,1), dydt])./2;
% calculate velocities and orientations
if ~isfield(simfile,'dT')
    simfile.dT = simfile.param.dT;
end
vx = dxdt./(simfile.dT*saveEvery);
vy = dydt./(simfile.dT*saveEvery);
speed = sqrt(vx(:,framesAnalyzed).^2 + vy(:,framesAnalyzed).^2);
N = simfile.N;
pairdist = NaN(N*(N-1)/2,numFrames);
dxcorr = NaN(size(pairdist));
vxcorr = NaN(size(pairdist));
distBins = 0:distBinWidth:maxDist;
gr = NaN(length(distBins) - 1,numFrames);
Area = simfile.L(1)*simfile.L(end); % should work for both scalar L and [Lx, Ly]
nNbrDist = NaN(N,numFrames);
for frameCtr = 1:numFrames
    frame = framesAnalyzed(frameCtr);
    if numel(trackedNodes)>1
        dxcorr(:,frameCtr) = vectorCrossCorrelation2D(ox(:,frame),oy(:,frame),true,true); % directional correlation
    end
    vxcorr(:,frameCtr) = vectorCrossCorrelation2D(vx(:,frame),vy(:,frame),true,false); % velocity correlation
    if strcmp(simfile.param.bc,'periodic')
        pairdist(:,frameCtr) = computeDistancesWithPeriodicBoundary([x(:,frame) y(:,frame)],simfile.L);
    else
        pairdist(:,frameCtr) = pdist([x(:,frame) y(:,frame)]); % distance between all pairs, in micrometer
    end
    gr(:,frameCtr) = histcounts(pairdist(:,frameCtr),distBins,'Normalization','count'); % radial distribution function
    gr(:,frameCtr) = gr(:,frameCtr)'.*Area./(2*pi*distBins(2:end)*distBinWidth)...
        /(N*(N-1)/2); % normalisation by number of pairs, not double-counting
    distanceMatrix = squareform(pairdist(:,frameCtr)); % distance of every worm to every other
    nNbrDist(:,frameCtr) = min(distanceMatrix + max(max(distanceMatrix))*eye(size(distanceMatrix)));
end
% bin distance data
[nNbrDistcounts,nNbrDistBins,nNbrDistbinIdx]  = histcounts(nNbrDist(:),...
    'BinWidth',distBinWidth,'BinLimits',[0 maxDist]);
[pairDistcounts,pairDistBins,pairDistbinIdx]  = histcounts(pairdist(:),...
    'BinWidth',distBinWidth,'BinLimits',[0 maxDist]);
% convert bin edges to centres (for plotting)
nNbrDistBins = nNbrDistBins(1:end-1) + diff(nNbrDistBins)/2;
pairDistBins = pairDistBins(1:end-1) + diff(pairDistBins)/2;
% ignore larger distance values and bins with only one element, as this will cause bootsci to fault
nNdistkeepIdcs = nNbrDistbinIdx>0&ismember(nNbrDistbinIdx,find(nNbrDistcounts>1));
nNbrDistBins = nNbrDistBins(nNbrDistcounts>1);
speed = speed(nNdistkeepIdcs);
nNbrDistbinIdx = nNbrDistbinIdx(nNdistkeepIdcs);
pdistkeepIdcs = pairDistbinIdx>0&ismember(pairDistbinIdx,find(pairDistcounts>1));
pairDistBins = pairDistBins(pairDistcounts>1);
dxcorr = dxcorr(pdistkeepIdcs);
vxcorr = vxcorr(pdistkeepIdcs);
pairDistbinIdx = pairDistbinIdx(pdistkeepIdcs);

[s_med,s_ci] = grpstats(speed(:),nNbrDistbinIdx,{@median,bootserr});
[corr_v_med,corr_v_ci] = grpstats(vxcorr(:),pairDistbinIdx,{@median,bootserr});
if numel(trackedNodes)>1
    [corr_o_med,corr_o_ci] = grpstats(dxcorr(:),pairDistbinIdx,{@median,bootserr});
else
    corr_o_med = NaN(size(corr_v_med));
    corr_o_ci = NaN(size(corr_v_ci));
end
end

function w = correctForPeriodicBoundary(v,L)
w = v;
w(v<=-L/2) = v(v<=-L/2) + L;
w(v>=L/2) = v(v>=L/2) - L;
end