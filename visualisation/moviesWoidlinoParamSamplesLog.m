% make movies from simulation results...

clear
close all

samplesToPlot = [1444];
nSamples = numel(samplesToPlot);
filepath = '/data1/linus/paramSamplesLog/woidlinos/';
addpath('../')

for sampleCtr = 1:nSamples
    thisSampleNum = samplesToPlot(sampleCtr);
    thisfile = dir([filepath '*v0_0.33_*sample_' num2str(thisSampleNum) '.mat']);
    if exist([filepath thisfile.name],'file')...
            &&~exist(['../movies/woidlinoMovies/paramSampleLogMovies/' strrep(thisfile.name,'.mat','.mp4')],'file')
        out = load([filepath thisfile.name]);
        animateWoidTrajectories(out.xyarray,['../movies/woidlinoMovies/paramSampleLogMovies/' strrep(thisfile.name,'.mat','.mp4')],...
            out.L,0.035,[],[0, 0]);
    elseif ~exist([filepath thisfile.name],'file')
        disp(['no results for ' thisfile.name])
    end
end
