function [ F_contact ] = resolveContacts(forceArray,distanceMatrixFull,distanceMatrix, objInd, nodeInd, cutoff)
% to resolve contact forces between overlapping nodes
% inputs:
% forceArray is N by M by ndim matrix of forces acting on every node
% distance matrix should have dimensions of N objects by M nodes by dimension
% objIdx is the scalar index of the agent upon which the force is to be
% calculated
% cutoff is the scalar radius below which volume exclusion is enforced
% outpus:
% F_contact will be column vector of dimension ndim

N = size(distanceMatrixFull,1);
M = size(distanceMatrixFull,2);
ndim = size(distanceMatrixFull,3);
collisionNbrs = distanceMatrix<=cutoff; % check distance too all other nodes of all other objects
collisionNbrs(objInd,max(nodeInd-1,1):min(nodeInd+1,M)) = false; % no contact force with self or adjacent nodes
% Nc = nnz(collisionNbrs);
if any(collisionNbrs(:))
    % find unit vectors pointing from neighbours to node
    if N>1
        e_nN = [distanceMatrixFull(collisionNbrs(:)) distanceMatrixFull(find(collisionNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
            ./repmat(distanceMatrix(collisionNbrs(:)),1,ndim); % normalise for distance
    elseif N==1 % special case due to how indexed entries return (annoying)
        e_nN = [distanceMatrixFull(collisionNbrs(:)) distanceMatrixFull(find(collisionNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
            ./repmat(distanceMatrix(collisionNbrs(:))',1,ndim); % normalise for distance
    end
    F_neighbours = diag([forceArray(collisionNbrs(:)) forceArray(find(collisionNbrs(:)) + N*M)]... % force acting on nodes in contact
        *e_nN'); %projected onto connecting vector: |fij*eij|
    % only resolve contact forces if force on neighbour is pointing towards
    % node in question
    F_neighbours(F_neighbours<0) = 0;
    F_contact = sum(repmat(F_neighbours,1,ndim).*e_nN,1); % contact force, summed over neighbours
else
    F_contact = [0; 0];
end
end

