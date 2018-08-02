function [] = plotActivityMap(filename,folder,nSnapShots)
% make a heatmap of worm positions from any simulation results...

% frameFraction is the % of frames from the end that get plotted
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',20,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1,...
    'Renderer','opengl');

addpath('../')
newFolder = '../figures/woids/activityMaps/';

activityFig = figure;
hold on
plotColors = parula(nSnapShots);
if exist([folder filename],'file')
    out = load([folder filename]);
    maxFrame = size(out.xyarray,4);
    framesSampled = round(linspace(1,maxFrame,nSnapShots+1));
    framesSampled = framesSampled(2:end); % removes first frame (= initial condition)
    for frameCtr = 1:nSnapShots
        thisFrame = framesSampled(frameCtr);
        positions2plot = out.xyarray(:,:,:,thisFrame);
        plotWoidTrajectoriesSingleFrame(positions2plot,out.L,0.035,...
            [plotColors(frameCtr,:) 0.5],false);
    end
    %% add colorbar
    hc = colorbar;
    hc.TickLabels = num2str(round(hc.Ticks'*60*out.T/3600)); % approx scale in minutes
    hc.Label.String = 'T (mins)';
    %% export figure
    figname = [newFolder strrep(filename,'.mat','.eps')];
    activityFig.PaperUnits = 'centimeters';
    exportfig(activityFig,figname, exportOptions)
    system(['epstopdf ' figname]);
    system(['rm ' figname]);
elseif ~exist([folder filename],'file')
    disp(['no results for ' filename])
end
