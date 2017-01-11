%% plot lag (or phase difference) of undulations along worm from single worm data
close all
clear

% find all files in the directory tree
files = rdir('../wormtracking/eigenworms/singleWorm/npr-1/**/*.mat');
numWorms = size(files,1);
framerate = 30;


velFig = figure; hold on
accFig = figure; hold on
for wormCtr = 1:numWorms
    % load single worm data
    load(files(wormCtr).name)
    [tangentAngles, ~] = makeAngleArrayV(worm.posture.skeleton.x',worm.posture.skeleton.y');
    % normalise for amplitude - divide by half the range
    tangentAngles = tangentAngles./repmat(range(tangentAngles)/2,size(tangentAngles,1),1);
    curvature = diff(tangentAngles,1,2);
    % estimate change and acceleration in angle over time at 30Hz
    dKdt = framerate*diff(curvature);
    d2Kdt2 = framerate*diff(dKdt);
    % calculate and plot group statistics
    Kbins = quant(curvature(1:end-1,:),0.01);
    [v, sigv] = grpstats(dKdt(:),Kbins(:),{@nanmean,@nanstd});
    boundedline(unique(Kbins(~isnan(Kbins(:)))),v,[sigv, sigv],...
        'alpha','transparency',1/numWorms,'nan','fill',velFig.Children)
    Kbins = quant(curvature(2:end-1,:),0.01);
    [a, siga] = grpstats(d2Kdt2(:),Kbins(:),{@mean,@std});
    boundedline(unique(Kbins(~isnan(Kbins(:)))),a,[siga, siga],...
        'alpha','transparency',1/numWorms,'nan','fill',accFig.Children)
end

%% edit plots
velFig.Children.XLabel.String = 'curvature (d\phi/ds)';
velFig.Children.YLabel.String = 'tangent angle velocity (d\phi/dt)';
velFig.Children.XLim = [-0.5 0.5];
velFig.Children.YLim = [-20 20];
% guide to the eye
x = -0.5:0.1:0.5; plot(velFig.Children,x,-20*x,'k--')

accFig.Children.XLabel.String = 'curvature (d\phi/ds)';
accFig.Children.YLabel.String = 'tangent angle acceleration (d^2\phi/dt^2)';
accFig.Children.XLim = [-0.5 0.5];
accFig.Children.YLim = [-1000 1000];
% guide to the eye
x = -0.5:0.1:0.5; plot(accFig.Children,x,-1200*x,'k--')
%% export figure
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',14,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

set(velFig,'PaperUnits','centimeters')
filename = ['parameterisationPlots/undulationStiffnessVel'];
exportfig(velFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

set(accFig,'PaperUnits','centimeters')
filename = ['parameterisationPlots/undulationStiffnessAcc'];
exportfig(accFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);