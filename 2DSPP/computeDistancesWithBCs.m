function distanceMatrix = computeDistancesWithBCs(positionArray,L,bc)
% find distances between all pairs of objects
ndim = size(positionArray,2);
N = size(positionArray,1);
distanceMatrix = NaN(N,N,ndim);
for dimCtr = 1:ndim
    distanceMatrix(:,:,dimCtr) = positionArray(:,dimCtr*ones(1,N)) - positionArray(:,dimCtr*ones(1,N))';
    % reset some distances if boundaries are periodic
    if (~iscell(bc)&&strcmp(bc,'periodic'))||(iscell(bc)&&strcmp(bc{dimCtr},'periodic'))
        if numel(L)==ndim % vector domain size [L_x L_y] ie non-square domain
            [overIndxRow, overIndxCol] = find(distanceMatrix(:,:,dimCtr)>=L(dimCtr)/2);
            distanceMatrix(overIndxRow, overIndxCol,dimCtr) = distanceMatrix(overIndxRow, overIndxCol,dimCtr) - L(dimCtr);
            [underIndxRow, underIndxCol] = find(distanceMatrix(:,:,dimCtr)<=-L(dimCtr)/2);
            distanceMatrix(underIndxRow, underIndxCol,dimCtr) = distanceMatrix(underIndxRow, underIndxCol,dimCtr) + L(dimCtr);
        else
            [overIndxRow, overIndxCol] = find(distanceMatrix(:,:,dimCtr)>=L/2);
            distanceMatrix(overIndxRow, overIndxCol,dimCtr) = distanceMatrix(overIndxRow, overIndxCol,dimCtr) - L;
            [underIndxRow, underIndxCol] = find(distanceMatrix(:,:,dimCtr)>=L/2);
            distanceMatrix(underIndxRow, underIndxCol,dimCtr) = distanceMatrix(underIndxRow, underIndxCol,dimCtr) + L;
        end
    end
end
end