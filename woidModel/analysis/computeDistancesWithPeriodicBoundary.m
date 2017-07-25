function distances = computeDistancesWithPeriodicBoundary(positionArray,L)
% find distances between all pairs of objects, correcting for periodic
% boundaries
% positionArray should have form [x, y], where x,y are column vectors

% issues / todo:

ndim = size(positionArray,2);
N = size(positionArray,1);
distanceMatrix = NaN(N,N,ndim);
for dimCtr = 1:ndim
    distanceMatrix(:,:,dimCtr) = bsxfun(@minus,positionArray(:,dimCtr),positionArray(:,dimCtr)');
    % reset some distances if boundaries are periodic
    if numel(L)==ndim % vector domain size [L_x L_y] ie non-square domain
        Ldim = L(dimCtr);
    else
        Ldim = L;
    end
    overIndcs = find(distanceMatrix(:,:,dimCtr)>=Ldim/2) + N^2*(dimCtr - 1);
    distanceMatrix(overIndcs) = distanceMatrix(overIndcs) - Ldim;
    underIndcs = find(distanceMatrix(:,:,dimCtr)<=-Ldim/2) + N^2*(dimCtr - 1);
    distanceMatrix(underIndcs) = distanceMatrix(underIndcs) + Ldim;
end
try
    distances = squareform(sqrt(sum(distanceMatrix.^2,3)));
catch
    1;
end
end