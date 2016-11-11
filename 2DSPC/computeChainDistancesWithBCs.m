function distanceMatrix = computeChainDistancesWithBCs(positionArray,L,bc)
% find distances between all pairs of objects

% issues / todo:
%   - very slow for periodic boundary conditions
ndim = size(positionArray,3);
N = size(positionArray,1);
M = size(positionArray,2);
distanceMatrix = NaN(N*M,N*M,ndim);
positionArrayStacked = reshape(positionArray,N*M,ndim);
for dimCtr = 1:ndim
    distanceMatrix(:,:,dimCtr) = positionArrayStacked(:,dimCtr*ones(1,N*M))...
        - positionArrayStacked(:,dimCtr*ones(1,N*M))';
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
distanceMatrix = reshape(distanceMatrix,N,M,N,M,ndim);
end