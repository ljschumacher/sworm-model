% make movies from simulation results...

clear
close all

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1,...
    'Renderer','painters');

samplesToPlot = [29339];
nSamples = numel(samplesToPlot);
filepath = '~/Dropbox/projects/collectiveBehaviour/sworm-model/woidModel/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r1/N2/';
addpath('../')

for sampleCtr = 1:nSamples
    snapshotFig = figure;
    thisSampleNum = samplesToPlot(sampleCtr);
    thisfile = dir([filepath '*v0_0.14_*sample_' num2str(thisSampleNum) '.mat']);
    if exist([filepath thisfile.name],'file')...
            &&~exist(['../figures/woidlinos/snapshots/' strrep(thisfile.name,'.mat','.pdf')],'file')
        out = load([filepath thisfile.name]);
        time2plot = round(size(out.xyarray,4)*(0.9 + 0.1*rand()));
        positions2plot = out.xyarray(:,:,:,time2plot);
        plotWoidTrajectoriesSingleFrame(positions2plot,out.L,0.035,[],true)
        figname = ['../figures/woidlinos/snapshots/' strrep(thisfile.name,'.mat','.eps')];
        %% export figure
        snapshotFig.PaperUnits = 'centimeters';
        exportfig(snapshotFig,figname, exportOptions)
        system(['epstopdf ' figname]);
        system(['rm ' figname]);
    elseif ~exist([filepath thisfile.name],'file')
        disp(['no results for ' thisfile.name])
    end
end
