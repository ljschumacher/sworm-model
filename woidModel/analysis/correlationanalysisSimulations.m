function [s_med,s_ci, corr_o_med,corr_o_ci, corr_v_med,corr_v_ci, gr,distBins] = ...
    correlationanalysisSimulations(simfile,trackedNodes,distBinWidth,framesAnalyzed)
% calculate directional & velocity correlation, and
% radial distribution functions

% issues/to-do:
% - for periodic bounadry conditions, if the tracked nodes overlap a
% boundary, it's not clear how to best find the centroid

% define functions for grpstats
mad1 = @(x) mad(x,1); % median absolute deviation
% alternatively could use boxplot-style confidence intervals on the mean,
% which are 1.57*iqr/sqrt(n) - unclear how justified this is
iqrci = @(x) 1.57*iqr(x)/sqrt(numel(x));
% or one could use a bootstrapped confidence interval
bootserr = @(x) bootci(1e2,{@nanmedian,x},'alpha',0.05,'Options',struct('UseParallel',true));

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
numFrames = numel(framesAnalyzed);
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
dxds = diff(simfile.xyarray(:,fliplr(trackedNodes),1,:),1,2);
dyds = diff(simfile.xyarray(:,fliplr(trackedNodes),2,:),1,2);
if strcmp(simfile.param.bc,'periodic')
    dxdt = correctForPeriodicBoundary(dxdt,simfile.L(1));
    dydt = correctForPeriodicBoundary(dydt,simfile.L(2));
    dxds = correctForPeriodicBoundary(dxds,simfile.L(1));
    dyds = correctForPeriodicBoundary(dyds,simfile.L(2));
end
% central difference (should be equibalent to matlab's gradient function
dxdt = ([dxdt, dxdt(:,end)] + [dxdt(:,1), dxdt])./2;
dydt = ([dydt, dydt(:,end)] + [dydt(:,1), dydt])./2;
% calculate velocities and orientations
if ~isfield(simfile,'dT')
    simfile.dT = simfile.param.dT;
end
vx = dxdt./(simfile.dT*simfile.saveevery);
vy = dydt./(simfile.dT*simfile.saveevery);
ox = squeeze(mean(dxds,2));
oy = squeeze(mean(dyds,2));
speed = sqrt(vx(:,framesAnalyzed).^2 + vy(:,framesAnalyzed).^2);
N = simfile.N;
pairdist = NaN(N*(N-1)/2,numFrames);
dxcorr = NaN(size(pairdist));
vxcorr = NaN(size(pairdist));
distBins = 0:distBinWidth:min(simfile.L)/2;
gr = NaN(length(distBins) - 1,numFrames);
nearestNbrDist = NaN(N,numFrames);
if numel(simfile.L)==1
    A = pi*simfile.L^2;
else
    A = simfile.L(1).*simfile.L(2);
end
for frameCtr = 1:numFrames
    frame = framesAnalyzed(frameCtr);
    dxcorr(:,frameCtr) = vectorCrossCorrelation2D(ox(:,frame),oy(:,frame),true); % directional correlation
    vxcorr(:,frameCtr) = vectorCrossCorrelation2D(vx(:,frame),vy(:,frame),true); % velocity correlation
    if strcmp(simfile.param.bc,'periodic')
        pairdist(:,frameCtr) = computeDistancesWithPeriodicBoundary([x(:,frame) y(:,frame)],simfile.L);
    else
        pairdist(:,frameCtr) = pdist([x(:,frame) y(:,frame)]); % distance between all pairs, in micrometer
    end
    gr(:,frameCtr) = histcounts(pairdist(:,frameCtr),distBins,'Normalization','count'); % radial distribution function
    gr(:,frameCtr) = gr(:,frameCtr)'.*A./(2*pi*distBins(2:end)*distBinWidth)...
        /N/(N-1); % normalisation
    distanceMatrix = squareform(pairdist(:,frameCtr)); % distance of every worm to every other
    nearestNbrDist(:,frameCtr) = min(distanceMatrix + max(max(distanceMatrix))*eye(size(distanceMatrix)));
end

[s_med,s_ci] = grpstats(speed(:),quant(nearestNbrDist(:),distBinWidth),...
    {@median,bootserr});
[corr_o_med,corr_o_ci] = grpstats(dxcorr(:),quant(pairdist(:),distBinWidth),...
    {@median,bootserr});
[corr_v_med,corr_v_ci] = grpstats(vxcorr(:),quant(pairdist(:),distBinWidth),...
    {@median,bootserr});
end

function w = correctForPeriodicBoundary(v,L)
w = v;
w(v<=-L/2) = v(v<=-L/2) + L;
w(v>=L/2) = v(v>=L/2) - L;
end