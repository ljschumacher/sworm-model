function distanceMatrix = computeDistancesWithBCs(positionArray,L,bc)
% find distances between all pairs of objects
ndim = size(positionArray,2);
N = size(positionArray,1);
distanceMatrix = NaN(N,N,ndim);
for dimCtr = 1:ndim
    distanceMatrix(:,:,dimCtr) = positionArray(:,dimCtr*ones(1,N)) - positionArray(:,dimCtr*ones(1,N))';
    % reset some distances if boundaries are periodic
    if (~iscell(bc)&&strcmp(bc,'periodic'))||(iscell(bc)&&strcmp(bc{dimCtr},'periodic'))
        if numel(L)==2 % vector domain size [L_x L_y] ie non-square domain
            [mirrorIndxRow, mirrorIndxCol] = find(abs(distanceMatrix(:,:,dimCtr))>=L(dimCtr)/2);
            distanceMatrix(mirrorIndxRow, mirrorIndxCol,dimCtr) = L(dimCtr) - ...
                distanceMatrix(mirrorIndxRow, mirrorIndxCol,dimCtr);
        else
            [mirrorIndxRow, mirrorIndxCol] = find(abs(distanceMatrix(:,:,dimCtr))>=L/2);
            distanceMatrix(mirrorIndxRow, mirrorIndxCol,dimCtr) = L - ...
                distanceMatrix(mirrorIndxRow, mirrorIndxCol,dimCtr);
        end
    end
end
end