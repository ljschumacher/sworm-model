% make movies from simulation results...

clear
close all

samplesToPlot = [17585 13108 14020];
nSamples = numel(samplesToPlot);
filepath = '/data2/linus/paramSamples/woidlinos/woidM18paramD2/';
addpath('../')

for sampleCtr = 1:nSamples
    thisSampleNum = samplesToPlot(sampleCtr);
    thisfile = dir([filepath '*v0_0.33_*sample_' num2str(thisSampleNum) '.mat']);
    if exist([filepath thisfile.name],'file')...
            &&~exist(['../movies/woidlinoMovies/paramSampleMovies/' strrep(thisfile.name,'.mat','.mp4')],'file')
        out = load([filepath thisfile.name]);
        animateWoidTrajectories(out.xyarray,['../movies/woidlinoMovies/paramSampleMovies/' strrep(thisfile.name,'.mat','.mp4')],...
            out.L,0.035,[],[0, 0]);
    elseif ~exist([filepath thisfile.name],'file')
        disp(['no results for ' thisfile.name])
    end
end
