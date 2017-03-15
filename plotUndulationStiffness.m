%% plot lag (or phase difference) of undulations along worm from single worm data
close all
clear

for strain = {'N2','npr-1','CB4856'}
    % find all files in the directory tree
    files = rdir(['../wormtracking/eigenworms/singleWorm/' strain{:} '/**/*.mat']);
    numWorms = min(size(files,1),100);
    framerate = 30;
    
    velFig = figure; hold on
    accFig = figure; hold on
    for wormCtr = 1:numWorms
        % load single worm data - somehow the filenames can have funny
        % characters in them from rdir which we need to remove
        load(strrep(files(wormCtr).name,'._',''))
        [segmentAngles, ~] = makeAngleArrayV(worm.posture.skeleton.x',worm.posture.skeleton.y');
        % normalise for amplitude - divide by half the range
        segmentAngles = segmentAngles./repmat(range(segmentAngles)/2,size(segmentAngles,1),1);
        deltaAngles = diff(segmentAngles,1,2);
        % estimate change and acceleration in angle over time at 30Hz
        dKdt = framerate*diff(deltaAngles);
        dK2dt = framerate*diff(diff(deltaAngles,1,2));
        % calculate and plot group statistics
        Kbins = quant(deltaAngles(1:end-1,:),0.01);
        [v, sigv] = grpstats(dKdt(:),Kbins(:),{@nanmean,@nanstd});
        boundedline(unique(Kbins(~isnan(Kbins(:)))),v,[sigv, sigv],...
            'alpha','transparency',1/numWorms,'nan','fill',velFig.Children)
        Kbins = quant(deltaAngles(1:end-1,1:end-1),0.02);
        [a, siga] = grpstats(dK2dt(:),Kbins(:),{@mean,@std});
        boundedline(unique(Kbins(~isnan(Kbins(:)))),a,[siga, siga],...
            'alpha','transparency',1/numWorms,'nan','fill',accFig.Children)
    end
    
    %% edit plots
    velFig.Children.XLabel.String = '\Delta\theta';
    velFig.Children.YLabel.String = 'intersegment angle change (d\theta/dt)';
    velFig.Children.XLim = [-0.5 0.5];
    velFig.Children.YLim = [-20 20];
    % guide to the eye
    x = -0.5:0.1:0.5; plot(velFig.Children,x,-20*x,'k--')
    title(velFig.Children,strain{:},'FontWeight','normal')
    
    accFig.Children.XLabel.String = '\Delta\theta';
    accFig.Children.YLabel.String = 'intersegment angle gradient change (d/dt(d\theta/ds))';
    accFig.Children.XLim = [-0.5 0.5];
    accFig.Children.YLim = [-20 20];
    % guide to the eye
    plot(accFig.Children,x,33*x,'k--')
    title(accFig.Children,strain{:},'FontWeight','normal')
    %% export figure
    exportOptions = struct('Format','eps2',...
        'Color','rgb',...
        'Width',14,...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',10,...
        'LineWidth',2);
    
    set(velFig,'PaperUnits','centimeters')
    filename = ['parameterisationPlots/undulationStiffnessVel' strain{:}];
    exportfig(velFig,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    system(['rm ' filename '.eps']);
    
    set(accFig,'PaperUnits','centimeters')
    filename = ['parameterisationPlots/undulationCurvatureVel' strain{:}];
    exportfig(accFig,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    system(['rm ' filename '.eps']);
end
tilefigs