function [] = snapshotsWoidlinoSimulation(filename,folder,frameFraction)
% make snapshot from any simulation results...

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1,...
    'Renderer','painters');

addpath('../')
newFolder = '../figures/woidlinos/snapshots/';

snapshotFig = figure;
if exist([folder filename],'file')
    out = load([folder filename]);
    maxFrame = size(out.xyarray,4);
    thisFrame = round(maxFrame*frameFraction);
    positions2plot = out.xyarray(:,:,:,thisFrame);
    plotWoidTrajectoriesSingleFrame(positions2plot,out.L,0.035,[],true)
    figname = [newFolder strrep(filename,'.mat','.eps')];
    %% export figure
    snapshotFig.PaperUnits = 'centimeters';
    exportfig(snapshotFig,figname, exportOptions)
    system(['epstopdf ' figname]);
    system(['rm ' figname]);
elseif ~exist([folder filename],'file')
    disp(['no results for ' filename])
end
