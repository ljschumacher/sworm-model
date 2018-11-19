% make movies from simulation results...

clear
close all

samplesToPlot = [5315];
nSamples = numel(samplesToPlot);
filepath = '~/Dropbox/projects/collectiveBehaviour/sworm-model/woidModel/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r2/npr_1/';
addpath('../')
moviepath = '../movies/woidlinoMovies/paramSampleMovies/r2/';
for sampleCtr = 1:nSamples
    thisSampleNum = samplesToPlot(sampleCtr);
    thisfile = dir([filepath '*v0_0.33_*sample_' num2str(thisSampleNum) '.mat']);
    if exist([filepath thisfile.name],'file')...
            &&~exist([moviepath strrep(thisfile.name,'.mat','.mp4')],'file')
        out = load([filepath thisfile.name]);
        animateWoidTrajectories(out.xyarray,[moviepath strrep(thisfile.name,'.mat','.mp4')],...
            out.L,0.035,[],[0, 0]);
    elseif ~exist([filepath thisfile.name],'file')
        disp(['no results for ' thisfile.name])
    end
end
