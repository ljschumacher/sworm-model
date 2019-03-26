% plot summary statistics of simulation and reference data


addpath('../../')
addpath('../')

clear
close all

% set pdf export options
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',6,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',9,...
    'LineWidth',1);

% set plotting options
plotColors = lines(2);
plotbins = (0.1:0.1:2) - 0.1/2;
strainCtr = 1; % which experiments to plot, 1 for npr-1, 2 for N2
ylabels = {'pcf','hcd'};

% set simulation file
filepath = '~/Dropbox/projects/collectiveBehaviour/sworm-model/results/woidlinos/paramSamples/PRW_4D_taxis_weighted_additive_r2/postiPredictiveCheck/';
figurepath = '../figures/woidlinos/sumstats/r2/';
filenames = rdir([filepath '*M45*run1.mat']);
nFiles = numel(filenames);
numReps = 3;

% set experimental summary statistics
sumstat_filename = '../../wormtracking/trackingAnalysis/inference/sumstats_expmnt.mat';
load(sumstat_filename);
% disregard tails of hierarchical clustering distribution
exp_ss_array{1,3} = exp_ss_array{1,3}(:,1:12);
exp_ss_array{2,3} = exp_ss_array{1,3}(:,1:12);

model_score = NaN(nFiles,numReps);

% loop over file to plot
for fileCtr = 1:nFiles
    fileName = strrep(filenames(fileCtr).name,filepath,'');
    shortFileName = strrep(strrep(strrep(strrep(strrep(fileName,'wlM18_N_40_L_7.5_v0_0.33_vs_0.018_',''),...
        'k_theta_0_slow_stochastic_bynode_dwell_0.0036_1.1_',''),...
        'revdensity_',''),'weighted_additive_',''),'_run1.mat','');
    if exist([filepath fileName],'file')
        for repCtr = 1:numReps
            thisFileName = strrep(fileName,'run1',['run' num2str(repCtr)]);
            out = load([filepath thisFileName]);
            %% analyse simulation data
            fraction_to_sample = min(out.param.dT*out.param.saveEvery/3,1); % specifiy fraction of frames to sample
            sim_sumstats{1}(repCtr,:) = inf_pcf(out.xyarray,'simulation',fraction_to_sample);
            sim_sumstats{2}(repCtr,:) = inf_hierarchicalclustering(out.xyarray,'simulation',fraction_to_sample);
            [sim_sumstats{3}(repCtr,:), sim_sumstats{4}(repCtr,:)] ...
           = inf_positionalmoments(out.xyarray, 'simulation', fraction_to_sample);
        end
        %% compute model score (distance used in inference)
        model_score(fileCtr,:) = f_model_score(exp_ss_array, sim_sumstats, [1 1 0 0]);
        disp(['score for ' shortFileName ' is ' num2str(mean(model_score(fileCtr,:),2),2)...
            '±' num2str(std(model_score(fileCtr,:),0,2)/sqrt(numReps),2)])
        %% plot summary stats
        for statCtr = 1:2
            sumstatFig = figure;
            nBins = size(exp_ss_array{strainCtr,statCtr+1},2);
            % plot experimental reference
            %             errorbar(plotbins(1:nBins),mean(exp_ss_array{strainCtr,statCtr+1}),...
            %                 std(exp_ss_array{strainCtr,statCtr+1}),':','LineWidth',1,'Color',plotColors(strainCtr,:))
            plot(plotbins(1:nBins),mean(exp_ss_array{strainCtr,statCtr+1}),...
                ':','LineWidth',1,'Color',plotColors(strainCtr,:))
            hold on
            sumstatFig.Children.YScale = 'log';
            % plot simulation result
            %             semilogy(plotbins(1:nBins),nanmean(sim_sumstats{statCtr},1),'LineWidth',1,'Color',plotColors(strainCtr,:))
            errorbar(plotbins(1:nBins),nanmean(sim_sumstats{statCtr},1),...
                nanstd(sim_sumstats{statCtr},0,1),'LineWidth',1,'Color',plotColors(strainCtr,:))
            if statCtr==2
                ylim([5e-3 1])
            end
            if statCtr==1
                ylim([0.5 100])
            end
            box on
%             grid on
            xlabel('r (mm)'), ylabel(ylabels{statCtr})
            figname = [figurepath 'S' num2str(statCtr) '_' shortFileName '.eps'];
            %% export figure
            sumstatFig.PaperUnits = 'centimeters';
            exportfig(sumstatFig,figname, exportOptions)
            system(['epstopdf ' figname]);
            system(['rm ' figname]);
        end
    elseif ~exist([filepath thisFileName],'file')
        disp(['no results for ' thisFileName])
    end
end