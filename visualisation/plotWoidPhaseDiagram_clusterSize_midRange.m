% plot woidlino phase portrait
% shows the end-points of of woidlet simulations for a 2D parameter sweep
close all
clear

addpath('../')
addpath('../analysis/')

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',17,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',1,...
    'Renderer','opengl');

M = 36;
L = [7.5 7.5];
N = 40;
plotColor = [0.25, 0.25, 0.25];
rc = 0.035;
ri = 5*rc;

revRatesClusterEdge = 0:5;

speed = [0.33];
slowspeed = 0.018;
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = 0:0.05:0.25;
dkdN_roam_values = 0:0.05:0.5;
% angleNoise = 1;
k_theta = 2;
% f_hapt = 0.5;

secondVariables = dkdN_dwell_values;
nrevRates = numel(revRatesClusterEdge);
ndwellVals = numel(dkdN_dwell_values);
nroamVals = numel(dkdN_roam_values);
aspectRatio = nrevRates/(ndwellVals);

BiggestComponentSize = NaN(nrevRates,ndwellVals,nroamVals);
rclust = 3*rc;
burnIn = 0.5;
%%
phasePortraitFig = figure;
% biggestComponentSizesAll = []; %%%%%%%%%
for revRateCtr = 1:nrevRates
    revRateClusterEdge = revRatesClusterEdge(revRateCtr);
    for ddwellCtr = 1:ndwellVals
        dkdN_dwell = dkdN_dwell_values(ddwellCtr);
        for droamCtr = 1:nroamVals
            dkdN_undwell = dkdN_roam_values(droamCtr);
            filename = ['woids_N_' num2str(N) '_L_' num2str(L(1)) ...
                ...'_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
                '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
                '_dkdN_' num2str(dkdN_dwell) '_' num2str(dkdN_undwell)... %% also check the case of indem dwell/roam
                '_revRateClusterEdge_' num2str(revRateClusterEdge,2) ...
                ...'_haptotaxis_' num2str(f_hapt) ...
                '_ri' num2str(ri) ...
                '_run1.mat'];
            filepath = '../results/woids/midRange/';
            if exist([filepath filename],'file')
                load([filepath filename])
                % loop over frames
                lastFrame = size(xyarray,4);
                firstFrame = round(lastFrame*burnIn);
                framesSampled = firstFrame:3:lastFrame;
                numSamples = numel(framesSampled);
                biggestComponentSizes = NaN(numSamples,1);
                for sampleCtr = 1:numSamples
                    thisFrame = framesSampled(sampleCtr);
                    % compute size of biggest connected component
                    biggestComponentSizes(sampleCtr) = calculateBiggestComponent(xyarray(:,:,:,thisFrame),rclust);
                end
                %                 if any(biggestComponentSizes>=30)
                % %                 biggestComponentSizesAll = [biggestComponentSizesAll,biggestComponentSizes];  %%%%%%%%
                %                    figure
                %                 animateWoidTrajectories(xyarray,...
                %                     ['../movies/woidMovies/mappingMovies/' filename(1:end-4)],L,0.035)
                % close(gcf)
                %                 end
                BiggestComponentSize(revRateCtr,ddwellCtr,droamCtr) = mean(biggestComponentSizes);
            else
                warning([filename ' does not exist'])
            end
        end
    end
end
% make diagonal slice to get coordinates for subsequent plotting
dslice = surf(dkdN_dwell_values,revRatesClusterEdge,dkdN_dwell_values.*ones(size(revRatesClusterEdge))');
xd = dslice.XData;
yd = dslice.YData;
zd = dslice.ZData;
delete(dslice)

h = slice(dkdN_dwell_values,revRatesClusterEdge,dkdN_roam_values,BiggestComponentSize,[0],[1 0],[ 0]);
hold on
h2 = slice(dkdN_dwell_values,revRatesClusterEdge,dkdN_roam_values,BiggestComponentSize,...
    repmat(xd,nrevRates,1),repmat(yd',1,ndwellVals),zd);
h = [h; h2];
nSlices = length(h);
for slcCtr = 1:nSlices
    h(slcCtr).FaceColor = 'interp';
    h(slcCtr).EdgeColor = 'none';
end
load ~/Dropbox/Utilities/colormaps_ascii/diverging/cmap_BuRd.txt
maxSize = ceil(max(BiggestComponentSize(:)));
colormap(cmap_BuRd(1:round(256/maxSize):end,:))
view(150,50)
ylabel('revRate')
xlabel('dk_{dwell}/d\rho')
zlabel('dk_{roam}/d\rho')
% plot an isosurface of a certain cluster size
% s = isosurface(dkdN_dwell_values,revRatesClusterEdge,dkdN_roam_values,BiggestComponentSize,30);
% add colorbar
hc = colorbar;
hc.Label.String = 'size of biggest cluster';
caxis([0 maxSize])
%% export figure
phasePortraitFig.PaperUnits = 'centimeters';
filename = ['../figures/woids/woidPhaseDiagram_clusterSize_3D'...
    '_N_' num2str(N) '_M_' num2str(M) '_L_' num2str(L(1)) '_noVolExcl' ...
    ...'_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
    '_speed_' num2str(speed,'%1.0e') ...
    '_slowing_' slowingMode '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
    '_ri' num2str(ri) ...
    '.eps'];
exportfig(phasePortraitFig,filename, exportOptions)
system(['epstopdf ' filename]);
system(['rm ' filename]);