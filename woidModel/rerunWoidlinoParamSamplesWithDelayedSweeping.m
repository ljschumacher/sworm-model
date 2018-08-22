function [] = rerunWoidlinoParamSamplesWithDelayedSweeping(sampleCtr,r_feed)

% rerun parameter samples with sweeping
addpath('visualisation/')

filepath = '~/Dropbox/projects/collectiveBehaviour/sworm-model/woidModel/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r1/npr_1/';

thisfile = dir([filepath '*v0_0.33_*sample_' num2str(sampleCtr) '.mat']);
if exist([filepath thisfile.name],'file')
    newfilename =  [strrep(thisfile.name,'.mat','') '_delayedsweeping_feedrate_' num2str(r_feed)];
    if~exist(['movies/woidlinoMovies/paramSampleMovies/' ...
            newfilename '.mp4'],'file')
        load([filepath thisfile.name]);
        % set parameters for sweeping simulation
        param.r_feed = r_feed;
        param.k_unroam = 10;
        Textend = 1000;
        % extend simulation - this will overwrite xyarray
        [xyarray, currentState, food] = runWoids(Textend,N,M,L,param,'resumeState',currentState);
        % save results, from continuation point only
        save(['results/woidlinos/paramSamples/sweeping/' newfilename '.mat'],...
            'xyarray','T','N','M','L','param','currentState','food')
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