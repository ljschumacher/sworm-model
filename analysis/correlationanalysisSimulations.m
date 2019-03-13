function [corr_o_mean,corr_o_ci, corr_v_mean,corr_v_ci,corr_vn_mean,corr_vn_ci, gr,distBins,nbrDistBins,pairDistBins] = ...
    correlationanalysisSimulations(simfile,trackedNodes,distBinWidth,framesAnalyzed,maxDist)
% calculate directional & velocity correlation, and
% radial distribution functions

% issues/to-do:
% - for periodic boundary conditions, if the tracked nodes overlap a
% boundary, it's not clear how to best find the centroid
% - Area in g(r) needs to be adjusted if arena is circular

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
    % normalise orientation for "head size"
    headSize = sqrt(ox.^2 + oy.^2);
    ox = ox./headSize;
    oy = oy./headSize;
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
N = simfile.N;
pairdist = NaN(N*(N-1)/2,numFrames);
dirxcorr = NaN(size(pairdist));
velxcorr = NaN(size(pairdist));
velnbrcorr = NaN(N*(N-1),numFrames);
distBins = 0:distBinWidth:maxDist;
gr = NaN(length(distBins) - 1,numFrames);
Area = simfile.L(1)*simfile.L(end); % should work for both scalar L and [Lx, Ly]
nbrDist = NaN(N*(N-1),numFrames);
for frameCtr = 1:numFrames
    frame = framesAnalyzed(frameCtr);
    if numel(trackedNodes)>1
        dirxcorr(:,frameCtr) = vectorCrossCorrelation2D(ox(:,frame),oy(:,frame),false,false); % directional correlation
    end
    velxcorr(:,frameCtr) = vectorCrossCorrelation2D(vx(:,frame),vy(:,frame),true,false); % velocity correlation
    dx = x(:,frame) - x(:,frame)';
    dy = y(:,frame) - y(:,frame)';
    if strcmp(simfile.param.bc,'periodic')
        pairdist(:,frameCtr) = computeDistancesWithPeriodicBoundary([x(:,frame) y(:,frame)],simfile.L);
        dx = correctForPeriodicBoundary(dx,simfile.L(1));
        dy = correctForPeriodicBoundary(dy,simfile.L(2));
    else
        pairdist(:,frameCtr) = pdist([x(:,frame) y(:,frame)]); % distance between all pairs, in micrometer
    end
    gr(:,frameCtr) = histcounts(pairdist(:,frameCtr),distBins,'Normalization','count'); % radial distribution function
    gr(:,frameCtr) = gr(:,frameCtr)'.*Area./(pi*(distBins(2:end).^2 - (distBins(2:end) - distBinWidth).^2))...
        /(N*(N-1)/2); % normalisation by number of pairs, not double-counting
    
    dD = sqrt(dx.^2 + dy.^2); % this should be equivalent to squareform(pairdist(:,frameCtr)
    selfNbrIndcs = dD(:)==0;
    velnbrcorrThisFrame = vectorPairedCorrelation2D(vx(:,frame),vy(:,frame),dx./dD,dy./dD,true,false);
    % don't keep entries where self is neighbour
    velnbrcorr(:,frameCtr) = velnbrcorrThisFrame(~selfNbrIndcs);
    nbrDist(:,frameCtr) = dD(~selfNbrIndcs);  
end
% bin distance data
[nbrDistcounts,nbrDistBins,nbrDistbinIdx]  = histcounts(nbrDist(:),...
    'BinWidth',distBinWidth,'BinLimits',[0 maxDist]);
[pairDistcounts,pairDistBins,pairDistbinIdx]  = histcounts(pairdist(:),...
    'BinWidth',distBinWidth,'BinLimits',[0 maxDist]);
% convert bin edges to centres (for plotting)
nbrDistBins = nbrDistBins + mean(diff(nbrDistBins))/2;
pairDistBins = pairDistBins + mean(diff(pairDistBins))/2;
% ignore larger distance values (bin=0) %and bins with only one element
nbrdistkeepIdcs = nbrDistbinIdx>0;%&ismember(nbrDistbinIdx,find(nbrDistcounts>1));
% nbrDistBins = nbrDistBins(nbrDistcounts>1);
nbrDistbinIdx = nbrDistbinIdx(nbrdistkeepIdcs);
velnbrcorr = velnbrcorr(nbrdistkeepIdcs);
pdistkeepIdcs = pairDistbinIdx>0;%&ismember(pairDistbinIdx,find(pairDistcounts>1));
% pairDistBins = pairDistBins(pairDistcounts>1);
dirxcorr = dirxcorr(pdistkeepIdcs);
velxcorr = velxcorr(pdistkeepIdcs);
pairDistbinIdx = pairDistbinIdx(pdistkeepIdcs);

[corr_v_mean,corr_v_ci] = grpstats(velxcorr(:),pairDistbinIdx,{'mean','meanci'});
[corr_vn_mean,corr_vn_ci] = grpstats(velnbrcorr(:),nbrDistbinIdx,{'mean','meanci'});

if numel(trackedNodes)>1
    [corr_o_mean,corr_o_ci] = grpstats(dirxcorr(:),pairDistbinIdx,{'mean','meanci'});
else
    corr_o_mean = NaN(size(corr_v_mean));
    corr_o_ci = NaN(size(corr_v_ci));
end
end

function w = correctForPeriodicBoundary(v,L)
w = v;
w(v<=-L/2) = v(v<=-L/2) + L;
w(v>=L/2) = v(v>=L/2) - L;
end