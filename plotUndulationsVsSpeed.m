%% plot undulation amplitude vs speed from single worm data

% find all files in the directory tree
files = rdir('../wormtracking/eigenworms/singleWorm/npr-1/**/*.mat');
numWorms = size(files,1);

% prealloc cell arrays to store wormwise data
amplitudes = cell(numWorms,1);
speeds = cell(numWorms,1);

for wormCtr = 1:numWorms
    % load single worm data
    load(files(wormCtr).name)
    
    speeds{wormCtr} = worm.locomotion.velocity.midbody.speed;
    amplitudes{wormCtr} = worm.locomotion.bends.midbody.amplitude;
end

%% plot
histogram2(horzcat(speeds{:}),horzcat(amplitudes{:}),'DisplayStyle','tile',...
    'Normalization','probability','EdgeColor','none')
xlabel('midbody speed')
ylabel('midbody bend amplitude')

%% export figure
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',14,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);
set(gcf,'PaperUnits','centimeters')
filename = ['parameterisationPlots/undulationsVspeed'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);