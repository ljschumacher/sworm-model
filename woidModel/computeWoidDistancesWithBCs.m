function distanceMatrix = computeWoidDistancesWithBCs(positionArray,L,bc)
% find distances between all pairs of objects

% issues / todo:
%   - very slow for periodic boundary conditions
%   - does not yet work with circular boundary
ndim = size(positionArray,3);
N = size(positionArray,1);
M = size(positionArray,2);
distanceMatrix = NaN(N*M,N*M,ndim);
positionArrayStacked = reshape(positionArray,N*M,ndim);
for dimCtr = 1:ndim
    distanceMatrix(:,:,dimCtr) = bsxfun(@minus,positionArrayStacked(:,dimCtr),positionArrayStacked(:,dimCtr)');
    % reset some distances if boundaries are periodic
    if (~iscell(bc)&&strcmp(bc,'periodic'))||(iscell(bc)&&strcmp(bc{dimCtr},'periodic'))
        if numel(L)==ndim % vector domain size [L_x L_y] ie non-square domain
            Ldim = L(dimCtr);
        else
            Ldim = L;
        end
        overIndcs = find(distanceMatrix(:,:,dimCtr)>=Ldim/2) + N^2*M^2*(dimCtr - 1);
        distanceMatrix(overIndcs) = distanceMatrix(overIndcs) - Ldim;
        underIndcs = find(distanceMatrix(:,:,dimCtr)<=-Ldim/2) + N^2*M^2*(dimCtr - 1);
        distanceMatrix(underIndcs) = distanceMatrix(underIndcs) + Ldim;
    end
end
distanceMatrix = reshape(distanceMatrix,N,M,N,M,ndim);
end