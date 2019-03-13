function [] = plotClusterSpeed()
% plot the speed of the CoM over time, for all simulations that start wuth
% filestring (these should only differ in their feeding rate)

filestring = '/Users/linus/Dropbox/projects/collectiveBehaviour/sworm-model/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r2/postiPredictiveCheck/wlM18_N_40_L_7.5_v0_0.33_vs_0.018_angleNoise_0.05_k_theta_0_slow_stochastic_bynode_dwell_0.0036_1.1_dkdN_0.32_0.98_revdensity_drdN_0.26_haptotaxis_weighted_additive_0.0019_postiMean_run1_sweeping_';

filenames = rdir([filestring '*']);
nFiles = numel(filenames);
nReps = 5;
feedrates = NaN(nFiles,1);
dt_experiment = 30;
com_velocity = NaN(nFiles,nReps,113);
for fileCtr = 1:nFiles
    for repCtr = 1:nReps
        if repCtr == 1
            filename = filenames(fileCtr).name;
        else
            filename = strrep(filenames(fileCtr).name,'run1',['run' num2str(repCtr)]);
        end
        try
        load(filename)
        dt = param.saveEvery*param.dT;
        % discard burn-in and "stationary/joining phase"
        xyarray = xyarray(:,:,:,round(2000/dt):end);
        [com_x, com_y] = centerOfMassPeriodicBoundary(xyarray,L);
        com_dx = diff(com_x);
        com_dy = diff(com_y);
        % filter for jumps in displacement >L/2
        com_dx(com_dx>L(1)/2) = com_dx(com_dx>L(1)/2) - L(1);
        com_dx(com_dx<-L(1)/2) = com_dx(com_dx<-L(1)/2) + L(1);
        com_dy(com_dy>L(2)/2) = com_dy(com_dy>L(2)/2) - L(2);
        com_dy(com_dy<-L(2)/2) = com_dy(com_dy<-L(2)/2) + L(2);
        % use this to 'unwrap' the original coordinates
        new_x = com_x(1) + [0; cumsum(com_dx)];
        new_y = com_y(1) + [0; cumsum(com_dy)];
        % subsample to reflect experimental framerate
        subsampleFactor = round(dt_experiment/dt);
        new_x = new_x(1:subsampleFactor:end);
        new_y = new_y(1:subsampleFactor:end);
        % calculate new displacements
        new_dx = diff(new_x);
        new_dy = diff(new_y);
        
        com_velocity(fileCtr,repCtr,:) = sqrt(new_dx.^2 + new_dy.^2)...
            ./(dt*subsampleFactor)*1000*60; % convert from mm/s to microns/min
        catch
            warning(['no file found for ' filename])
        end
    end
    feedrates(fileCtr) = param.r_feed;
end
% combine results from all replicates
com_velocity = reshape(com_velocity,nFiles,[]);
% plot results
clustSpeedFig = figure;
errorbar(feedrates*100,median(com_velocity/100,2),mad(com_velocity/100,1,2),'LineWidth',2)
hold on
% plot experimental reference
plot([0 1],1.72*ones(1,2),'k--')
xlim([0 1])
ylim([0 5])
xlabel('feeding rate / lawn thickness (relative units)')
ylabel('cluster centroid speed (100\mum/min)')

%% export plot
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',12,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1);
clustSpeedFig.PaperUnits = 'centimeters';
filename = ['../figures/woidlinos/' strrep(filestring,[filenames(1).folder '/'],'') ...
    'clusterSpeed.eps'];
exportfig(clustSpeedFig,filename, exportOptions)
system(['epstopdf ' filename]);
system(['rm ' filename]);
end

function [com_x, com_y] = centerOfMassPeriodicBoundary(xyarray,L)
% this calculates the centre of mass for peridic boundaries - useful
% trick found on wikipedia
c_x = mean(mean(cos(xyarray(:,:,1,:)/L(1)*2*pi),2),1);
s_x = mean(mean(sin(xyarray(:,:,1)/L(1)*2*pi),2),1);
com_x = squeeze(L(2)/2/pi*(atan2(-s_x,-c_x) + pi) - L(1)/2);
c_y = mean(mean(cos(xyarray(:,:,2,:)/L(2)*2*pi),2),1);
s_y = mean(mean(sin(xyarray(:,:,2,:)/L(2)*2*pi),2),1);
com_y = squeeze(L(2)/2/pi*(atan2(-s_y,-c_y) + pi) - L(2)/2);
end