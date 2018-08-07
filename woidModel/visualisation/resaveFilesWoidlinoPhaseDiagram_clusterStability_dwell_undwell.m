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
    'Renderer','painters');

M = 18;
L = [3.6 3.6];
N = 40;
plotColor = [0.25, 0.25, 0.25];

revRatesClusterEdge = 0:1:5;

speed = [0.33];
% slowspeeds = fliplr([0.33, 0.05, 0.025, 0.0125]);
slowspeed = 0.018;
slowingMode = 'stochastic_bynode';
k_dwell = 0.0036;
k_undwell = 1.1;
dkdN_dwell_values = 0:0.1:1;
dkdN_roam_values = 0:0.2:2;
% angleNoise = 1;
k_theta = 2;
% f_hapt = 0.5;

secondVariables = dkdN_dwell_values;
nrevRates = numel(revRatesClusterEdge);
ndwellVals = numel(dkdN_dwell_values);
nroamVals = numel(dkdN_roam_values);
aspectRatio = nrevRates/(ndwellVals);

Rgyr = NaN(nrevRates,ndwellVals,nroamVals); % radius of gyration
ri = 3*0.035;
%%
for revRateCtr = 1:nrevRates
    revRateClusterEdge = revRatesClusterEdge(revRateCtr);
    for ddwellCtr = 1:ndwellVals
        dkdN_dwell = dkdN_dwell_values(ddwellCtr);
        for droamCtr = 1:nroamVals
            dkdN_undwell = dkdN_roam_values(droamCtr);
            oldFilename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
                ...'_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
                '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
                '_dkdN_' num2str(dkdN_dwell) '_' num2str(dkdN_undwell)...
                '_revRateClusterEdge_' num2str(revRateClusterEdge,1) ...
                ...'_haptotaxis_' num2str(f_hapt) ...
                '_clusteredStart' ...
                '_run1.mat'];
            filepath = '../results/woidlinos/';
            if exist([filepath oldFilename],'file')
                load([filepath oldFilename])
                if param.revRateClusterEdge~=revRateClusterEdge
                    disp(['resaving' oldFilename])
                    newFilename = ['wlM' num2str(M) '_N_' num2str(N) '_L_' num2str(L(1)) ...
                        ...'_angleNoise_' num2str(angleNoise) '_k_theta_' num2str(k_theta)...
                        '_v0_' num2str(speed,'%1.0e') '_vs_' num2str(slowspeed,'%1.0e') ...
                        '_' slowingMode 'SlowDown' '_dwell_' num2str(k_dwell) '_' num2str(k_undwell)...
                        '_dkdN_' num2str(dkdN_dwell) '_' num2str(dkdN_undwell)...
                        '_revRateClusterEdge_' num2str(param.revRateClusterEdge,2) ...
                        ...'_haptotaxis_' num2str(f_hapt) ...
                        '_clusteredStart' ...
                        '_run1'];
                    save([filepath newFilename '.mat'],'xyarray','T','N','M','L','param','currentState')
                    delete([filepath oldFilename])
                end
            end
        end
    end
end