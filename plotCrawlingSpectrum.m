%% plot the orthogonal-to-CoM velocities for different worm segments from single worm data

% issues / todo:
% - the N2 spectra average each other out
clear
close all
for strain = {'N2','npr-1','CB4856'}
    % find all files in the directory tree
    files = rdir(['../wormtracking/eigenworms/singleWorm/' strain{:} '/**/*.mat']);
    numWorms = min(size(files,1),100);
    numSegments = 49;
    
    % prealloc cell arrays to store wormwise data
    dN = cell(numWorms,1);
    dP = cell(numWorms,1);
    worm_freq = NaN(numWorms,1);
    for wormCtr = 1:numWorms
        % load single worm data - somehow the filenames can have funny
        % characters in them from rdir which we need to remove
        load(strrep(files(wormCtr).name,'._',''))
        
        % identify frames with missing data - find out error codes for missing
        % data in info.video.annotations
        nanframes = any(isnan(worm.posture.skeleton.x),1)|any(isnan(worm.posture.skeleton.x),1);
        numSamples = nnz(~nanframes);
        
        % convert positions into components of displacement
        dx = diff(worm.posture.skeleton.x,1,2);
        dy = diff(worm.posture.skeleton.y,1,2);
        
        % average along the worm body to get CoM motion
        dX = mean(dx); dY = mean(dy);
        
        % check if CoM speed is constant: plot(sqrt(dX.^2 + dY.^2))
        
        %project velocity components [dx, dy] to orthogonal of CoM, i.e. [-dY, dX] (rotated
        %90 degree counter clockwise)
        dn = (-dx.*dY(ones(numSegments,1),:) + dy.*dX(ones(numSegments,1),:)) ...% broadcasting using Tony's trick
            ./sqrt(dX(ones(numSegments,1),:).^2 + dY(ones(numSegments,1),:).^2); % normalised by CoM speed
        %project onto parallel to CoM
        dp = (dx.*dX(ones(numSegments,1),:) + dy.*dY(ones(numSegments,1),:)) ...% broadcasting using Tony's trick
            ./sqrt(dX(ones(numSegments,1),:).^2 + dY(ones(numSegments,1),:).^2); % normalised by CoM speed
        % normalise by speed
        dn0 = dn./sqrt(dp.^2 + dn.^2);
        dp0 = dp./sqrt(dp.^2 + dn.^2);
        % check forward, paused, reverse moving states
        fwdIdx = worm.locomotion.motion.mode==1;
        %         psdIdx = worm.locomotion.motion.mode==0;
        %         bwdIdx = worm.locomotion.motion.mode==-1;
        
        dN{wormCtr} = dn0(:,fwdIdx(1:end-1));
        dP{wormCtr} = dp0(:,fwdIdx(1:end-1));
        % %     plot(find(fwdIdx(1:end-1)),dn0(1,fwdIdx(1:end-1)))
        
        worm_freq(wormCtr) = nanmean(abs(worm.locomotion.bends.head.frequency(fwdIdx)));
    end
    
    %% plot angles of movement of all worms combined
    % histogram2(horzcat(dP{:}),horzcat(dN{:}),'DisplayStyle','tile')
    % hold on, plot(0,0,'r+') % mark origin
    
    % %plot the angles of velocity vector to CoM
    % histogram(atan2(horzcat(dN{:}),horzcat(dP{:})),'EdgeColor','none')
    
    % dNall = horzcat(dN{:});
    % dPall = horzcat(dP{:});
    % figure, hold on
    % for segCtr = 1:12:49
    %     histogram(atan2(dNall(segCtr,:),dPall(segCtr,:)),'DisplayStyle','stairs','Normalization','probability')
    % end
    
    %% plot the power spectrum of tangential velocity or angle time series
    % averaged over all worms
    segCtr = 1;
    Fs = 25; % sampling frequency
    % preallocate variables
    [~, ncols] = cellfun(@size,dN);
    Lmax = max(ncols);
    Lmax = Lmax + mod(Lmax,2); % make even
    PN_nan0 = NaN(numWorms,Lmax/2+1);
    PN_itp1 = NaN(numWorms,Lmax/2+1);
    PN_itp3 = NaN(numWorms,Lmax/2+1);
    f = Fs*(0:(Lmax/2))/Lmax; % vector of frequencies for plotting
    for wormCtr = 1:numWorms
        % select time series of chosen body segment
        dn = dN{wormCtr}(segCtr,:);
        L = numel(dn);
        if mod(L,2)
            dn(end+1) = 0;
            L = L+1;
        end
        t_idx = 1:L;
        nan_idx = isnan(dn);
        if any(~nan_idx)
            % variant a) set NaNs = 0
            dn_nan0 = dn; dn_nan0(nan_idx) = 0;
            % variant b) linearly interpolate NaNs
            dn_itp1 = dn; dn_itp1(nan_idx) = ...
                interp1(t_idx(~nan_idx),dn_itp1(~nan_idx),t_idx(nan_idx),'linear',0);
            % variant c) cubic interpolate NaNs
            dn_itp3 = dn; dn_itp3(nan_idx) = ...
                interp1(t_idx(~nan_idx),dn_itp3(~nan_idx),t_idx(nan_idx),'pchip',0);
            % compute fourier transforms and power spectra
            % a
            fn_nan0 = fft(dn_nan0);
            Pn_nan0 = abs(fn_nan0/L); % two-sided power spectrum
            Pn_nan0 = Pn_nan0(1:L/2+1); Pn_nan0(2:end-1) = 2*Pn_nan0(2:end-1);% single-sided power spectrum
            PN_nan0(wormCtr,1:(L/2+1)) = Pn_nan0;
            % b
            fn_itp1 = fft(dn_itp1);
            Pn_itp1 = abs(fn_itp1/L); % two-sided power spectrum
            Pn_itp1 = Pn_itp1(1:L/2+1); Pn_itp1(2:end-1) = 2*Pn_itp1(2:end-1);% single-sided power spectrum
            PN_itp1(wormCtr,1:(L/2+1)) = Pn_itp1;
            % c
            fn_itp3 = fft(dn_itp3);
            Pn_itp3 = abs(fn_itp3/L); % two-sided power spectrum
            Pn_itp3 = Pn_itp3(1:L/2+1); Pn_itp3(2:end-1) = 2*Pn_itp3(2:end-1);% single-sided power spectrum
            PN_itp3(wormCtr,1:(L/2+1)) = Pn_itp3;
        end
    end
    %% plot
    wsize = 5; % for moving average smoothing power spectrum
    figure
    plot(f,movmean(nanmean(PN_itp3),wsize))
    hold on
    plot(f,movmean(nanmean(PN_itp1),wsize))
    plot(f,movmean(nanmean(PN_nan0),wsize))
    % plot the frequency from database for comparison
    plot(mean(worm_freq).*[1 1],[0 0.05],'k--')
    legh = legend('cubic','linear','NaN=0');
    legh.Title.String = 'interpolation';
    legh.Title.FontWeight = 'normal';
    xlabel('frequency (Hz)')
    ylabel('amplitude^2')
    ylim([0 0.04])
    xlim([0 3])
    title(strain{:},'FontWeight','normal')
    %% export figure
    exportOptions = struct('Format','eps2',...
        'Color','rgb',...
        'Width',10,...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',10,...
        'LineWidth',2);
    set(gcf,'PaperUnits','centimeters')
    filename = ['parameterisationPlots/crawlingSpectrum' strain{:}];
    exportfig(gcf,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    system(['rm ' filename '.eps']);
end