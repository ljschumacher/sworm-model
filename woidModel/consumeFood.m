function foodGrid = consumeFood(foodGrid,xgrid,ygrid,r_feed,dT,positions)
% consume food from grid based on worm head positions
% issues/to-do:
% - currently only works for rectangular domains
headPositions_x = positions(:,1,1);
headPositions_y = positions(:,1,2);
N = size(positions,1);
% for each worm-head, find the closest grid point
[~, feedx] = min(abs(headPositions_x - xgrid),[],2);
[~, feedy] = min(abs(headPositions_y - ygrid),[],2);
% consume food from the grid
for objCtr = 1:N
foodGrid(feedx(objCtr),feedy(objCtr)) =...
    max(foodGrid(feedx(objCtr),feedy(objCtr)) - r_feed*dT,0);
end
end

