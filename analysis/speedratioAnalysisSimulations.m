function [speedFractions, densitybinedges] = ...
    speedratioAnalysisSimulations(simfile,trackedNodes,framesAnalyzed,nBins,speedThresh)
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

if nargin<4||isempty(framesAnalyzed)
    framesAnalyzed = 1:size(simfile.xyarray,4);
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
end
% central difference (should be equivalent to matlab's gradient function
dxdt = ([dxdt, dxdt(:,end)] + [dxdt(:,1), dxdt])./2;
dydt = ([dydt, dydt(:,end)] + [dydt(:,1), dydt])./2;
% calculate velocities and orientations
if ~isfield(simfile,'dT')
    simfile.dT = simfile.param.dT;
end
vx = dxdt./(simfile.dT*saveEvery)*1000;
vy = dydt./(simfile.dT*saveEvery)*1000;
N = simfile.N;
knn6density = NaN(N,numFrames);
speed = sqrt(vx(:,framesAnalyzed).^2 + vy(:,framesAnalyzed).^2);
for frameCtr = 1:numFrames
    frame = framesAnalyzed(frameCtr);
    dx = x(:,frame) - x(:,frame)';
    dy = y(:,frame) - y(:,frame)';
    if strcmp(simfile.param.bc,'periodic')
        dx = correctForPeriodicBoundary(dx,simfile.L(1));
        dy = correctForPeriodicBoundary(dy,simfile.L(2));
    end
    
    distanceMatrix = sqrt(dx.^2 + dy.^2); % this should be equivalent to squareform(pairdist(:,frameCtr)
    % sort distances
    distanceMatrix = sort(distanceMatrix,2);
    knn6density(:,frameCtr) = 6./(pi*distanceMatrix(:,6+1).^2); % shift index by 1 as closest neighbour is self
end
% 2d histcounts of absolute speed vs density
[~, densitybinedges, ~] = histcounts2(knn6density,speed,...
    'BinWidth',[0.5, 10],'XBinLimits',[0 10],'Normalization','Probability');

for binCtr = 1:nBins
    thisBinLogInd = knn6density>=densitybinedges(binCtr)&knn6density<densitybinedges(binCtr+1);
    [thisbincounts, thisbinedges] = histcounts(speed(thisBinLogInd),...
        'BinWidth',10,'Normalization','Probability');
    % calculate proportion in high and low speed states
    speedFractions(binCtr,1) = sum(thisbincounts(thisbinedges(2:end)<=speedThresh));
    speedFractions(binCtr,2) = sum(thisbincounts(thisbinedges(1:end-1)>speedThresh));
end

end

function w = correctForPeriodicBoundary(v,L)
w = v;
w(v<=-L/2) = v(v<=-L/2) + L;
w(v>=L/2) = v(v>=L/2) - L;
end