function model_score = f_model_score(exp_ss_array, sim_sumstats, weights, strainCtr)
% Compute the appropriate distances between each of the
% simulations and the experimental references


num_statistics = size(exp_ss_array,2)-1;

if nargin<3
    weights=[];
end
if nargin<4
    strainCtr = 1;
end

numSims = size(sim_sumstats{1},1);
expsim_dists = zeros(numSims, 1+num_statistics);

for statCtr = 1:num_statistics
    exp_data = exp_ss_array{strainCtr,1+statCtr};
    assert(~any(exp_data(:)==0),'zero-experimental data, cannot compute log-distances')
    for simCtr = 1:numSims
        if statCtr>2
            1;
        end
        sim_data = sim_sumstats{statCtr}(simCtr,:);
        dim_factor = 1./sqrt(size(exp_data,2)); % correction factor for higher dimensional summary statistics
        % Compute the distance between this simulation and the
        % reference - careful not to take log(0)
        if size(exp_data,2)>1
            sim_data = max(sim_data,5e-3); % to prevent log(0)
        end
        expsim_dists(simCtr,1+statCtr) = sum(vecnorm(...
            (log(exp_data) - log(sim_data))... % take scaled difference of all observed values of this summary stat and this simulated one
            ,2,2)... % take norm for each expmntl sample
            .*dim_factor... % correct for dim of summary stat
            );% sum this distance over expmntl samples
        
        % add the distance to the total from all summary statistics
        expsim_dists(simCtr,1) = expsim_dists(simCtr,1)...
            + weights(statCtr).*expsim_dists(simCtr,1+statCtr); % weight summary statistic
        if statCtr==num_statistics&&expsim_dists(simCtr,1)==0
            warning('zero distance btw expmnt and simulation');
        end
    end
end

model_score = expsim_dists(:,1);
end