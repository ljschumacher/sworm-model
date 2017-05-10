function [ F_contact ] = resolveContacts(forceArray,distanceMatrixFull,distanceMatrix, objInd, nodeInd, r_collision, sigma_LJ, r_LJcutoff, eps_LJ)
% to resolve contact forces between overlapping nodes
% inputs:
% forceArray is N by M by ndim matrix of forces acting on every node
% distance matrix should have dimensions of N objects by M nodes by dimension
% objIdx is the scalar index of the agent upon which the force is to be
% calculated
% cutoff is the scalar radius below which volume exclusion is enforced
% outpus:
% F_contact will be column vector of dimension ndim

% issues / to-dos:
% - should woids have contact forces with themselves?

N = size(distanceMatrixFull,1);
M = size(distanceMatrixFull,2);
ndim = size(distanceMatrixFull,3);
collisionNbrs = distanceMatrix<r_collision; % check distance to all other nodes of all other objects
collisionNbrs(objInd,:) = false; % no contact force with self or for adjacent nodes: max(nodeInd-1,1):min(nodeInd+1,M)
if nargin < 8
    eps_LJ = 0;
    if nargin < 7
        r_LJcutoff = 0;
    end
end
attractionNbrs = distanceMatrix<=r_LJcutoff;
attractionNbrs(collisionNbrs) = false;
attractionNbrs(objInd,:) = false;
% contact forces
if any(collisionNbrs(:))
    % find unit vectors pointing from neighbours to node
    e_nN = [distanceMatrixFull(collisionNbrs(:)) distanceMatrixFull(find(collisionNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
        ./distanceMatrix(collisionNbrs(:)); % normalise for distance
    F_neighbours = diag([forceArray(collisionNbrs(:)) forceArray(find(collisionNbrs(:)) + N*M)]... % force acting on nodes in contact
        *e_nN'); % force of neighbour projected onto connecting vector: |fij*eij|
    % only resolve contact forces if force on neighbour is pointing towards
    % node in question
    F_neighbours(F_neighbours<0) = 0;
    F_contact = sum(F_neighbours.*e_nN,1); % contact force, summed over neighbours
else
    F_contact = [0; 0];
end
% adhesion forces
if any(attractionNbrs(:))&&N>1
    e_nN = [distanceMatrixFull(attractionNbrs(:)) distanceMatrixFull(find(attractionNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
        ./distanceMatrix(attractionNbrs(:)); % normalise for distance
    f_LJ = 48*eps_LJ./distanceMatrix(attractionNbrs(:)).*((sigma_LJ./distanceMatrix(attractionNbrs(:))).^12 ...
        - 1/2*(sigma_LJ./distanceMatrix(attractionNbrs(:))).^6);
    F_LJ = sum(f_LJ.*e_nN,1); % adhesion force, summed over neighbours
    if ~any(collisionNbrs(:))
        F_contact = F_contact + F_LJ';
    else
        F_contact = F_contact + F_LJ;
    end
end

end

