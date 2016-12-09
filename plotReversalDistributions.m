%% plot distribution of reversal events from single worm data

% find all files in the directory tree
files = rdir('../wormtracking/eigenworms/singleWorm/N2/**/*.mat');
numWorms = size(files,1);

% prealloc cell arrays to store wormwise data
numRevs = NaN(numWorms,1);
numFrames = NaN(numWorms,1);
% % figure, hold on
for wormCtr = 1:numWorms
    % load single worm data
    load(files(wormCtr).name)
    
    % count number of reversal periods
    numRevs(wormCtr) = size(worm.locomotion.motion.backward.frames,2);
    numFrames(wormCtr) = size(worm.locomotion.motion.mode,2);
end

%% save variables
save('reversalEvents.mat','numRevs','numFrames','numWorms')
%% plot
histogram(numRevs./numFrames,'Normalization','probability',...
    'EdgeColor','none')
xlim([0 0.01])
xlabel('events')
ylabel('P')

%% export figure
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',14,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);
set(gcf,'PaperUnits','centimeters')
filename = ['parameterisationPlots/reversalDistribution'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);