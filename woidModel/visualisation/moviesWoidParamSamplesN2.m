% make movies from simulation results...

clear
close all

samplesToPlot = [1000 4689 74];
nSamples = numel(samplesToPlot);
% filepath = '/datafast/linus/paramSampleResults/woids/';
filepath = ['/datafast/linus/paramSampleResults/woids/'];
addpath('../')

for sampleCtr = 1:nSamples
    thisSampleNum = samplesToPlot(sampleCtr);
    thisfile = dir([filepath '*v0_0.14_*sample_' num2str(thisSampleNum) '.mat']);
    if exist([filepath thisfile.name],'file')...
            &&~exist(['../movies/woidMovies/paramSampleMovies/' strrep(thisfile.name,'.mat','.mp4')],'file')
        out = load([filepath thisfile.name]);
        animateWoidTrajectories(out.xyarray,['../movies/woidMovies/paramSampleMovies/' strrep(thisfile.name,'.mat','.mp4')],...
            out.L,0.035,[],[0, 0]);
    elseif ~exist([filepath thisfile.name],'file')
        disp(['no results for ' thisfile.name])
    end
end
