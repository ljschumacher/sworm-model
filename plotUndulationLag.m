%% plot lag (or phase difference) of undulations along worm from single worm data
close all
clear

% find all files in the directory tree
files = rdir('../wormtracking/eigenworms/singleWorm/npr-1/**/*.mat');
numWorms = size(files,1);

% prealloc cell arrays to store wormwise data
deltaPhase = cell(numWorms,1);
figure, hold on
for wormCtr = 1:numWorms
    % load single worm data
    load(files(wormCtr).name)
    [tangentAngles, ~] = makeAngleArrayV(worm.posture.skeleton.x',worm.posture.skeleton.y');
    % normalise for amplitude and get phase
    % dphase = 2asin(sqrt(2mean((deltaPhi/2/amplitude).^2)))
    deltaPhase{wormCtr} = 2*asin(sqrt(2*diff(...
        tangentAngles./repmat(range(tangentAngles)/2,size(tangentAngles,1),1)/2,...
        1,2).^2));
    plot(nanmean(deltaPhase{wormCtr}))
end

%% plot
xlabel('worm segment')
ylabel('estimated phase difference btw segments')

%% export figure
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',14,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);
set(gcf,'PaperUnits','centimeters')
filename = ['parameterisationPlots/phaseDifferenceAlongWorm'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);