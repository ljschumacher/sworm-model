function [] = rerunWoidlinoParamSamplesWithSweeping(sampleCtr,r_feed)

% rerun parameter samples with sweeping
addpath('visualisation/')

filepath = '~/Dropbox/projects/collectiveBehaviour/sworm-model/woidModel/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r1/npr_1/';

thisfile = dir([filepath '*v0_0.33_*sample_' num2str(sampleCtr) '.mat']);
if exist([filepath thisfile.name],'file')
    newfilename =  [strrep(thisfile.name,'.mat','') '_sweeping_feedrate_' num2str(r_feed)];
    if~exist(['movies/woidlinoMovies/paramSampleMovies/' ...
            newfilename '.mp4'],'file')
        load([filepath thisfile.name]);
        % set parameters for sweeping simulation
        param.r_feed = r_feed;
        param.k_unroam = 10;
        T = 4000;
        param.saveEvery = round(2/param.dT);
        % run simulation
        [xyarray, currentState, food] = runWoids(T,N,M,L,param);
        save(['results/woidlinos/paramSamples/sweeping/' newfilename '.mat'],'xyarray','T','N','M','L','param','currentState')
        % make movie
        animateWoidTrajectories(xyarray,['movies/woidlinoMovies/paramSampleMovies/' newfilename '.mp4'],...
            L,0.035,food);
    elseif exist(['movies/woidlinoMovies/paramSampleMovies/' ...
            newfilename '.mp4'],'file')
        disp(['sweeping movie already exists for ' thisfile.name])
    end
elseif ~exist([filepath thisfile.name],'file')
    disp(['no results for ' thisfile.name])
end

end