function [ F_contact ] = resolveContactsLoop(forceArray,distanceMatrixFull,distanceMatrix,...
    r_collision, sigma_LJ, r_LJcutoff, eps_LJ, LJnodes, LJmode)
% to resolve contact forces between overlapping nodes
% inputs:
% forceArray is N by M by ndim by N by M matrix of forces acting on every node
% distance matrix should have dimensions of N objects by M nodes by N by M
% objIdx is the scalar index of the agent upon which the force is to be
% calculated
% cutoff is the scalar radius below which volume exclusion is enforced
% outpus:
% F_contact will be column vector of dimension ndim

% issues / to-dos:
% - should woids have contact forces with themselves?

N = size(distanceMatrixFull,1);
M = size(distanceMatrixFull,2);
F_contact = NaN(N,M,2);
F_LJ = NaN(N,M,2);

for objCtr = 1:N
    for nodeCtr = 1:M
        thisDistanceMatrix = distanceMatrix(:,:,objCtr,nodeCtr);
        thisDistanceMatrixFull = distanceMatrixFull(:,:,:,objCtr,nodeCtr);
        collisionNbrs = thisDistanceMatrix<r_collision; % check distance to all other nodes of all other objects
        collisionNbrs(objCtr,:) = false; % no contact force with self or for adjacent nodes: max(nodeInd-1,1):min(nodeInd+1,M)
        % contact forces
        if any(collisionNbrs(:))
            % find unit vectors pointing from neighbours to node
            e_nN = bsxfun(@rdivide,[thisDistanceMatrixFull(collisionNbrs(:)) thisDistanceMatrixFull(find(collisionNbrs(:)) + N*M)],...
                thisDistanceMatrix(collisionNbrs(:))); % bsxfun has similar performace to implicit expansion (below) but is mex-file compatible
            %     e_nN = [distanceMatrixFull(collisionNbrs(:)) distanceMatrixFull(find(collisionNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
            %         ./distanceMatrix(collisionNbrs(:)); % normalise for distance
            F_neighbours = diag([forceArray(collisionNbrs(:)) forceArray(find(collisionNbrs(:)) + N*M)]... % force acting on nodes in contact
                *e_nN'); % force of neighbour projected onto connecting vector: |fij*eij|
            % only resolve contact forces if force on neighbour is pointing towards
            % node in question
            F_neighbours(F_neighbours<0) = 0;
            F_contact(objCtr,nodeCtr,:) = sum(bsxfun(@times,F_neighbours,e_nN),1); % contact force, summed over neighbours
        else
            F_contact(objCtr,nodeCtr,:) = [0; 0];
        end
        % adhesion forces
        if ismember(nodeCtr,LJnodes)
            lennardjonesNbrs = thisDistanceMatrix<=r_LJcutoff;
            lennardjonesNbrs(objCtr,:) = false;
            if any(lennardjonesNbrs(:))&&N>1
                e_nN = bsxfun(@rdivide,[thisDistanceMatrixFull(lennardjonesNbrs(:)) thisDistanceMatrixFull(find(lennardjonesNbrs(:)) + N*M)],...
                    thisDistanceMatrix(lennardjonesNbrs(:))); % bsxfun has similar performace to implicit expansion (below) but is mex-file compatible
                %   e_nN = [distanceMatrixFull(lennardjonesNbrs(:)) distanceMatrixFull(find(lennardjonesNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
                %         ./distanceMatrix(lennardjonesNbrs(:)); % normalise for distance
                switch LJmode
                    case 'hard'
                        asigr = (0.2599210499*sigma_LJ + thisDistanceMatrix(lennardjonesNbrs(:)));
                        f_LJ = 8*eps_LJ./asigr.*((sigma_LJ./asigr).^4 - 1/2*sigma_LJ./asigr);
                    case 'soft'
                        asigr = (2/3*sigma_LJ + thisDistanceMatrix(lennardjonesNbrs(:)));
                        f_LJ = 8*eps_LJ./asigr.*((sigma_LJ./asigr).^2 - 1/2*sigma_LJ./asigr);
                    otherwise
                        f_LJ = zeros(size(lennardjonesNbrs(:)));
                end
                F_LJ(objCtr,nodeCtr,:) = sum(bsxfun(@times,f_LJ,e_nN),1); % adhesion force, summed over neighbours
            else
                F_LJ(objCtr,nodeCtr,:) = [0; 0];
            end
        else
            F_LJ(objCtr,nodeCtr,:) = [0; 0];
        end
    end
end
F_contact = F_contact + F_LJ;
end