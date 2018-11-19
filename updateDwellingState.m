function dwellLogInd = updateDwellingState(dwellLogInd,k_dwell,k_undwell,dkdN_dwell,dkdN_undwell,num_nbr_per_node,dT)
% updates movement state of worm based on specified rates of
% entering/leaving the dwelling state

% fast worms start dwelling with rate k_dwell + dkdN_dwell*num_nbr_per_node
if any(~dwellLogInd)
    speeders = find(~dwellLogInd);
    dwellrate =  k_dwell + dkdN_dwell.*num_nbr_per_node;
    starters = logical(poissrnd(dwellrate(speeders)*dT,numel(speeders),1));
    dwellLogInd(speeders(starters)) = true;
end

% dwelling worms stop dwelling with rate k_undwell, decreasing
%   exponentially with more neighbours
% NOTE: thsi current implentation allows worms that just stopped to speed
% up again. you may not want this, but it's currently kept for backwards
% compatibility (consistency with previous results) 
if any(dwellLogInd)
    dwellers = find(dwellLogInd);
    undwellrate = k_undwell*exp(-dkdN_undwell.*num_nbr_per_node);
    stoppers = logical(poissrnd(undwellrate(dwellers)*dT,numel(dwellers),1));
    dwellLogInd(dwellers(stoppers)) = false;
end

end

