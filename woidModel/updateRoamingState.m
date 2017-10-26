function roamingLogInd = updateRoamingState(roamingLogInd,k_roam,k_unroam, dT)
% updates movement state of worm based on specified rates of
% entering/leaving the roaming state

roamers = find(roamingLogInd);
nonroamers = find(~roamingLogInd);

% non-roaming worms start roaming with rate k_roam
starters = logical(poissrnd(k_roam*dT,numel(nonroamers),1));
roamingLogInd(nonroamers(starters)) = true;

% roaming worms stop roaming with rate k_unroam
stoppers = logical(poissrnd(k_unroam*dT,numel(roamers),1));
roamingLogInd(roamers(stoppers)) = false;

end

