function [ F_contact ] = resolveContacts(forceArray,distanceMatrixFull,distanceMatrix,...
    objInd, nodeInd, r_collision, sigma_LJ, r_LJcutoff, eps_LJ, LJnodes)
% to resolve contact forces between overlapping nodes
% inputs:
% forceArray is N by M by ndim matrix of forces acting on every node
% distance matrix should have dimensions of N objects by M nodes
% objIdx is the scalar index of the agent upon which the force is to be
% calculated
% cutoff is the scalar radius below which volume exclusion is enforced
% outpus:
% F_contact will be column vector of dimension ndim

% issues / to-dos:
% - should woids have contact forces with themselves?

N = size(distanceMatrixFull,1);
M = size(distanceMatrixFull,2);
% ndim = size(distanceMatrixFull,3);
if nargin < 8
    eps_LJ = 0;
    if nargin < 7
        r_LJcutoff = 0;
    end
end
collisionNbrs = distanceMatrix<r_collision; % check distance to all other nodes of all other objects
collisionNbrs(objInd,:) = false; % no contact force with self or for adjacent nodes: max(nodeInd-1,1):min(nodeInd+1,M)
% contact forces
if any(collisionNbrs(:))
    % find unit vectors pointing from neighbours to node
    e_nN = bsxfun(@rdivide,[distanceMatrixFull(collisionNbrs(:)) distanceMatrixFull(find(collisionNbrs(:)) + N*M)],...
        distanceMatrix(collisionNbrs(:))); % bsxfun has similar performace to implicit expansion (below) but is mex-file compatible
    %     e_nN = [distanceMatrixFull(collisionNbrs(:)) distanceMatrixFull(find(collisionNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
    %         ./distanceMatrix(collisionNbrs(:)); % normalise for distance
    F_neighbours = diag([forceArray(collisionNbrs(:)) forceArray(find(collisionNbrs(:)) + N*M)]... % force acting on nodes in contact
        *e_nN'); % force of neighbour projected onto connecting vector: |fij*eij|
    % only resolve contact forces if force on neighbour is pointing towards
    % node in question
    F_neighbours(F_neighbours<0) = 0;
    F_contact = sum(bsxfun(@times,F_neighbours,e_nN),1); % contact force, summed over neighbours
else
    F_contact = [0; 0];
end
% adhesion forces
if ismember(nodeInd,LJnodes) % check if current node feels LJ force
    lennardjonesNbrs = distanceMatrix<=r_LJcutoff;
    lennardjonesNbrs(objInd,:) = false;
    if any(lennardjonesNbrs(:))&&N>1
        e_nN = bsxfun(@rdivide,[distanceMatrixFull(lennardjonesNbrs(:)) distanceMatrixFull(find(lennardjonesNbrs(:)) + N*M)],...
            distanceMatrix(lennardjonesNbrs(:))); % bsxfun has similar performace to implicit expansion (below) but is mex-file compatible
        %   e_nN = [distanceMatrixFull(lennardjonesNbrs(:)) distanceMatrixFull(find(lennardjonesNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
        %         ./distanceMatrix(lennardjonesNbrs(:)); % normalise for distance
        f_LJ = 48*eps_LJ./distanceMatrix(lennardjonesNbrs(:)).*((sigma_LJ./distanceMatrix(lennardjonesNbrs(:))).^12 ...
            - 1/2*(sigma_LJ./distanceMatrix(lennardjonesNbrs(:))).^6);
        F_LJ = sum(bsxfun(@times,f_LJ,e_nN),1); % adhesion force, summed over neighbours
        if ~any(collisionNbrs(:))
            F_contact = F_contact + F_LJ';
        else
            F_contact = F_contact + F_LJ;
        end
    end
end
end

