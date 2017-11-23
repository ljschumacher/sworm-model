function roamingLogInd = updateRoamingState(roamingLogInd,k_roam,k_unroam,dT,...
    foodGrid,xgrid,ygrid,r_feed,positions)
% updates movement state of worm based on specified rates of
% entering/leaving the roaming state, as well as presence of food

roamers = find(roamingLogInd);
nonroamers = find(~roamingLogInd);

% non-roaming worms start roaming with rate k_roam
starters = logical(poissrnd(k_roam*dT,numel(nonroamers),1));
roamingLogInd(nonroamers(starters)) = true;

% worms that have run out of food start roaming
if r_feed>0
    headPositions_x = positions(:,1,1);
    headPositions_y = positions(:,1,2);
    N = size(positions,1);
    % for each worm-head, find the closest grid point
    [~, feedx] = min(abs(headPositions_x - xgrid),[],2);
    [~, feedy] = min(abs(headPositions_y - ygrid),[],2);
    % check food on the grid
    for objCtr = 1:N
        if foodGrid(feedx(objCtr),feedy(objCtr)) <= eps
            roamingLogInd(objCtr) = true;
        end
    end
end
% roaming worms stop roaming with rate k_unroam
stoppers = logical(poissrnd(k_unroam*dT,numel(roamers),1));
roamingLogInd(roamers(stoppers)) = false;

end

