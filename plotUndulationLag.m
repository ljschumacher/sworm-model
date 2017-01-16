%% plot lag (or phase difference) of undulations along worm from single worm data
close all
clear

for strain = {'npr-1','CB4856','N2'} 
    % find all files in the directory tree
    files = rdir(['../wormtracking/eigenworms/singleWorm/' strain{:} '/**/*.mat']);
    numWorms = min(size(files,1),100);
    meanDeltaPhase = NaN(numWorms,1);
    % prealloc cell arrays to store wormwise data
    figure, hold on
    for wormCtr = 1:numWorms
        % load single worm data - somehow the filenames can have funny
        % characters in them from rdir which we need to remove
        load(strrep(files(wormCtr).name,'._',''))
        [tangentAngles, ~] = makeAngleArrayV(worm.posture.skeleton.x',worm.posture.skeleton.y');
        % normalise for amplitude and get phase
        % dphase = 2asin(sqrt(2mean((deltaPhi/amplitude/2).^2)))
        deltaPhase = 2*asin(sqrt(2*diff(...
            tangentAngles./repmat(range(tangentAngles)/2,size(tangentAngles,1),1)/2,... % amplitude = range(tangentAngles)/2
            1,2).^2));
        plot(nanmean(deltaPhase))
        meanDeltaPhase(wormCtr) = mean(nanmean(deltaPhase));
    end
    
    %% plot
    plot([0, 50],mean(meanDeltaPhase)*[1, 1],'k--','LineWidth',2)
    xlabel('worm segment')
    ylabel('estimated phase difference btw segments')
    title(strain{:},'FontWeight','normal')
    ylim([0 0.3])
    %% export figure
    exportOptions = struct('Format','eps2',...
        'Color','rgb',...
        'Width',14,...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',10,...
        'LineWidth',2);
    set(gcf,'PaperUnits','centimeters')
    filename = ['parameterisationPlots/phaseDifferenceAlongWorm' strain{:}];
    exportfig(gcf,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    system(['rm ' filename '.eps']); 
end