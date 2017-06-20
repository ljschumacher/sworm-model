function [s_med,s_mad, corr_o_med,corr_o_mad, corr_v_med,corr_v_mad, gr,distBins] = ...
    correlationanalysisSimulations(simfile,trackedNodes,distBinWidth,framesAnalyzed)
% calculate directional & velocity correlation, and
% radial distribution functions

% issues/to-do:

% define functions for grpstats
mad1 = @(x) mad(x,1); % median absolute deviation
% alternatively could use boxplot-style confidence intervals on the mean,
% which are 1.57*iqr/sqrt(n)
M = size(simfile.xyarray,2);
if nargin<2||isempty(trackedNodes)
    trackedNodes = 1:M;
end
if nargin<3||isempty(distBinWidth)
    distBinWidth = 0.01; % in units of mm
end
if nargin<4||isempty(framesAnalyzed)
    framesAnalyzed = 1:size(simfile.xyarray,4);
end
numFrames = numel(framesAnalyzed);
%% calculate stats
x = squeeze(mean(simfile.xyarray(:,trackedNodes,1,:),2)); % centroid of tracked obj
y = squeeze(mean(simfile.xyarray(:,trackedNodes,2,:),2)); % centroid of tracked obj
vx = gradient(x,simfile.dT*simfile.saveevery);
vy = gradient(y,simfile.dT*simfile.saveevery);
ox = squeeze(mean(diff(simfile.xyarray(:,fliplr(trackedNodes),1,:),1,2),2));
oy = squeeze(mean(diff(simfile.xyarray(:,fliplr(trackedNodes),2,:),1,2),2));
speed = sqrt(vx(:,framesAnalyzed).^2 + vy(:,framesAnalyzed).^2);
N = simfile.N;
pairdist = NaN(N*(N-1)/2,numFrames);
dxcorr = NaN(size(pairdist));
vxcorr = NaN(size(pairdist));
distBins = 0:distBinWidth:min(simfile.L);
gr = NaN(length(distBins) - 1,numFrames);
nearestNbrDist = NaN(N,numFrames);
for frameCtr = 1:numFrames
    frame = framesAnalyzed(frameCtr);
    dxcorr(:,frameCtr) = vectorCrossCorrelation2D(ox(:,frame),oy(:,frame),true); % directional correlation
    vxcorr(:,frameCtr) = vectorCrossCorrelation2D(vx(:,frame),vy(:,frame),true); % velocity correlation
    pairdist(:,frameCtr) = pdist([x(:,frame) y(:,frame)]); % distance between all pairs, in micrometer
    gr(:,frameCtr) = histcounts(pairdist(:,frameCtr),distBins,'Normalization','count'); % radial distribution function
    gr(:,frameCtr) = gr(:,frameCtr)'.*simfile.L^2./(2*distBins(2:end)*distBinWidth)...
        /N/(N-1); % normalisation
    distanceMatrix = squareform(pairdist(:,frameCtr)); % distance of every worm to every other
    nearestNbrDist(:,frameCtr) = min(distanceMatrix + max(max(distanceMatrix))*eye(size(distanceMatrix)));
end
[s_med,s_mad] = grpstats(speed(:),quant(nearestNbrDist(:),distBinWidth),...
    {@median,mad1});
[corr_o_med,corr_o_mad] = grpstats(dxcorr(:),quant(pairdist(:),distBinWidth),...
    {@median,mad1});
[corr_v_med,corr_v_mad] = grpstats(vxcorr(:),quant(pairdist(:),distBinWidth),...
    {@median,mad1});
