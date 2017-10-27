function dwellLogInd = updateDwellingState(dwellLogInd,k_dwell,k_undwell,dT)
% updates movement state of worm based on specified rates of
% entering/leaving the dwelling state

dwellers = find(dwellLogInd);
speeders = find(~dwellLogInd);

% fast worms start dwelling with rate k_dwell
starters = logical(poissrnd(k_dwell*dT,numel(speeders),1));
dwellLogInd(speeders(starters)) = true;

% dwelling worms stop dwelling with rate k_undwell
stoppers = logical(poissrnd(k_undwell*dT,numel(dwellers),1));
dwellLogInd(dwellers(stoppers)) = false;

end

