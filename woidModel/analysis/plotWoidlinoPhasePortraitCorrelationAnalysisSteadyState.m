function [] = plotWoidlinoPhasePortraitCorrelationAnalysisSteadyState
% plot phase portrait of g(r) vs time to show quasi steady state

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

N = 40;
M = 18;
L = 7.5;
numRepeats = 1;
trackedNodes = 1:max(round(M*0.16),1);
maxDist = 1.2;
distBinwidth = 0.1;

% -- slow-down parameters --
vs = 0.018;% vs: speed when slowed down (default v0/3)
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
% -- reversal parameters --
reversalMode = 'density';
% -- Lennard-Jones parameters --
r_LJcutoff = -1;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
sigma_LJ = 0;  % particle size for Lennard-Jones force
eps_LJ = 0;
% -- undulation parameters --
k_theta = 0;
theta_0 = 0;
omega_m = 0;
deltaPhase = 0;
angleNoise = 0.05;
% -- haptotaxis --
% f_hapt = 0.5;
% -- speed and time-step --
v0 = 0.33; % npr1 0.33; N2 0.14

filepath = '../results/woidlinos/paramSamples/';
files = dir([filepath '*.mat']);
nPlots = 100;

pcfFig = figure;
for plotCtr = 1:nPlots
    thisFile = load([filepath files(plotCtr).name]);
    maxNumFrames = size(thisFile.xyarray,4);
    burnIn = round(1000./thisFile.T*maxNumFrames); % used for visualizing cut-off
    numFrames =  0.1*maxNumFrames; % sample some percentage of frames
    framesAnalyzed = round(linspace(1,maxNumFrames,numFrames)); % regularly sample frames without replacement
    
    %% calculate stats
    [~,~, ~,~, ~,~,...
        gr,~,~,~] = ...
        correlationanalysisSimulations(thisFile,trackedNodes,distBinwidth,framesAnalyzed,maxDist);
    % get maximum of gr
    [grmax, ~] = max(smoothdata(gr,2,'movmean',3));
    % plot lines for this file
    % speed v distance
    sumStat = smoothdata(grmax,'movmean',3);
    plot(framesAnalyzed,sumStat./mean(sumStat(round(burnIn./maxNumFrames*numel(sumStat)):end)),...
        'Color',[0 0 0 0.25])
    if plotCtr==1
        hold on
        plot(burnIn*[1 1],[0 max(grmax)],'k--')
    end
end

ylim([0 2])
ylabel('summart statistic (relative to mean after burn-in')
xlim([0 maxNumFrames])
xlabel('timesteps')
%% export figures
% radial distribution / pair correlation
pcfFig.PaperUnits = 'centimeters';
fignameprefix = ['figures/diagnostics/grmaxOverTime'];
fignamesuffix = ['_wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
    '_v0_' num2str(v0) '_vs_' num2str(vs) ...
    '_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
    '_slow_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
    '_rev' reversalMode  ...
    '_' num2str(nPlots) 'samples.eps'];
filename = [fignameprefix 'Radialdistribution' fignamesuffix];
exportfig(pcfFig,filename, exportOptions)
system(['epstopdf ' filename]);
system(['rm ' filename]);
end