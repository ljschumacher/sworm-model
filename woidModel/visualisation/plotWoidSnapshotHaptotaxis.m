% plot woid phase portrait
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

radius = 0.035;
plotColor = [0.25, 0.25, 0.25];
N = 40;
Lval = 7.5;
revRatesClusterEdge = 0:5;
speeds = [0.33];
% slowspeeds = fliplr([0.33, 0.1, 0.05, 0.025, 0.0125, 0.005]);
% slowspeeds = fliplr([0.33, 0.025, 0.0125, 0.005, 0.001]);
slowspeeds = [0.018];
slowingMode = 'stochastic_bynode';
eps_LJ = 1e-3;
f_hapt = 0.1;

k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = 0:0.2:1;
secondVariables = dkdN_dwell_values;
nrevRates = numel(revRatesClusterEdge);
ndwellVals = numel(dkdN_dwell_values);
aspectRatio = nrevRates/(ndwellVals + 1/3);
numRepeats = 1;
maxMeanClusterSize = 10;
burnIn = 0.9;
rclust = 3*radius;
for repCtr=1:numRepeats
    for speed = speeds
        snapshotFig = figure;
        for slowspeed = slowspeeds
            for dkdN_dwell = dkdN_dwell_values
                for revRateClusterEdge = revRatesClusterEdge
                    filename = ['woids_N_' num2str(N) '_L_' num2str(Lval) ...
                        '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                        '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell) ...
                        '_dkdN_' num2str(dkdN_dwell) '_' num2str(dkdN_dwell)...
                        '_revRateClusterEdge_' num2str(revRateClusterEdge,2)...
                        ...'_LJsoft' num2str(eps_LJ) ...
                        '_haptotaxis_' num2str(f_hapt) ...
                        '_run' num2str(repCtr) '.mat'];
                    filepath = '../results/woids/mapping/';
                    if exist([filepath filename],'file')
                        load([filepath filename],'xyarray')
                        % calculate average cluster size
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
                        if mean(biggestComponentSizes)>=maxMeanClusterSize
                            bestFile = filename;
                            maxMeanClusterSize = mean(biggestComponentSizes);
                            disp(['best sim so far is ' bestFile ' with average cluster size of '...
                                num2str(maxMeanClusterSize)])
                        end
                    else
                        disp(['No file ' filename])
                    end
                end
            end
        end
        %% plot selected simulation outcome
        snapshotsWoidSimulation(bestFile,filepath,0.95)
    end
end