% Function to calculate the pair correlation function over many sampled
% frames, and return a discretised array of the resulting g(r) distribution
function pcf_mean = inf_pcf(data, format, fraction_to_sample)

% Create bins of a given width, to store the data in. The bin width
% controls the width of the rings drawn from each reference particle during
% the computation of g(r). It also controls the bins into which the g(r)
% distribution is discretised.
bin_width = 0.1;
L = 7.5;
bins = 0:bin_width:1.2;

if nargin<3
    fraction_to_sample = 0.1; % specify the proportion of frames to be sampled
end

if strcmp(format,'simulation') || strcmp(format,'complexsim')||strcmp(format,'simulation-test')
    burn_in = 0.5; % specifies how much to ignore at the start of the simulation
    
    % Get the dimensions of the dataframe
    dims = size(data);
    if strcmp(format,'simulation')
        trackedNodes = 1:3;% only track nodes equivalent to the head
    elseif strcmp(format,'complexsim')
        trackedNodes = 1:6;
    elseif strcmp(format,'simulation-test')
        trackedNodes = 1;
    end
    
    % Get information from the dimensions of the input data
    num_worms = dims(1);
    final_frame = dims(4);
    
    % Sample fraction of the frames in the video
    num_samples = max(3,ceil(final_frame * (1 - burn_in) * fraction_to_sample));
    sampled_frames = randi([round(burn_in*final_frame) final_frame],1,num_samples); % could also sample without replacement, or regularly, but it doesn't seem to make a difference
    for sampleCtr = 1:num_samples
        
        % Access the data for the tracked worm node(s)
        thisFrameData = data(:,round(mean(trackedNodes)),:,sampled_frames(sampleCtr));
        
        % Initialise empty matrices to store location coordinates
        coords = zeros(num_worms,2);
        
        coords(:,1) = thisFrameData(:,:,1);
        coords(:,2) = thisFrameData(:,:,2);
        
        % Calculate pairwise distances with custom distance function
        % 'periodiceucdists' to take into account the horizontal and
        % vertical periodicity of the simulations.
        pair_dist = pdist(coords, @periodiceucdist);
        
        % Get the histogram counts of the pair_dist data using the bins
        pcf_raw = histcounts(pair_dist,bins,'Normalization','count');
        
        % Radial distribution function
        % Normalization step
        pcf_normalised = pcf_raw.*L^2./(pi*(bins(2:end).^2 - (bins(2:end) - bin_width).^2)*num_worms*(num_worms-1)/2); % normalisation by number of pairs, not double-counting
        
        % Store the gr information for each of the sampled timepoints
        if sampleCtr == 1
            pcf_store = zeros(num_samples,length(pcf_normalised));
        end
        pcf_store(sampleCtr,:) = pcf_normalised;
    end
    
    % Compute the average g(r) over the sampled timepoints
    pcf_mean = mean(pcf_store);
    
elseif strcmp(format,'experiment')
    % Analagous code for obtaining the same gr output from the
    % experimental .hdf5 data structure
    
    % Make pixel xy coordinates informative by converting to mm
    pix2mm = 0.1/19.5;
    
    % Randomly sample fraction of the frames in the video
    frames = data{3};
    num_samples = max(3,ceil(length(unique(frames)) * fraction_to_sample));
    frames_sampled = randi([min(frames),max(frames)], 1, num_samples);
    
    for sampleCtr = 1:num_samples
        thisFrame = frames_sampled(sampleCtr);
        
        thisFrame_ind = find(frames==thisFrame);
        
        while length(thisFrame_ind) < 2 % resample if less than two worms in frame
            thisFrame = randi([min(frames),max(frames)],1);
            frames_sampled(sampleCtr) = thisFrame;
            thisFrame_ind = find(frames==thisFrame);
        end
        
        num_worms = length(thisFrame_ind);
        coords = zeros(num_worms,2);
        
        coords(:,1) = data{1}(thisFrame_ind).*pix2mm;
        coords(:,2) = data{2}(thisFrame_ind).*pix2mm;
        
        % Obtain the pairwise distances
        pair_dist = pdist(coords);
        
        % Get the histogram counts of the pair_dist data using the bins
        pcf_raw = histcounts(pair_dist,bins,'Normalization','count');
        
        % Radial distribution function
        % Normalization step
        pcf_normalised = pcf_raw.*(pi*(8.5/2).^2)...
            ./(pi*(bins(2:end).^2 - (bins(2:end) - bin_width).^2)*num_worms*(num_worms-1)/2);

        % Store the gr information for each of the sampled timepoints
        if sampleCtr == 1
            pcf_store = zeros(num_samples,length(pcf_normalised));
        end
        pcf_store(sampleCtr,:) = pcf_normalised;
    end
    
    % Compute the average g(r) over the sampled timepoints
    pcf_mean = mean(pcf_store);
end
end