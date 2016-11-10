function [ F_excl ] = exclusionForce(distanceMatrixFull, objIdx, cutoff)
% calculates repulsive exclusion force between agents separated by less
% that cutoff distance
% inputs:
% distance matrix should have dimensions of N by N agents by dimension
% objIdx is the scalar index of the agent upon which the force is to be
% calculated
% cutoff is the scalar radius below which volume exclusion is enforced
% outpus:
% F_excl will be column vector of same dimension

ndim = size(distanceMatrixFull,3);
distanceMatrix = sqrt(sum(distanceMatrixFull.^2,3)); % reduce to scalar
collisionNbrs = distanceMatrix(:,objIdx)<=cutoff;
collisionNbrs(objIdx) = false; % no self-repulsion
Nc = nnz(collisionNbrs);
    F_excl = sum(... % sum over all collision neighbours
        reshape(distanceMatrixFull(objIdx,collisionNbrs,:),Nc,ndim)... %direction FROM neighbours TO object (use reshape rather than squeeze for case when Nc = 1)
        ./distanceMatrix(collisionNbrs,objIdx*ones(1,ndim))... % normalise for distance
        ,1)'; % made into column vector

end

