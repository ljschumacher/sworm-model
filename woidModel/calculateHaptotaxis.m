function [Fh] = calculateHaptotaxis(distanceMatrixFull,distanceMatrix,objInd,ri,r_overlap,f_hapt)
% calculates the haptotactic contribution to the motile force, i.e. the
% bias of movement towards the mean direction of Nbrs within an interaction
% radius
N = size(distanceMatrixFull,1);
M = size(distanceMatrixFull,2);
% check which nodes are within interaction radius
haptoNbrs = distanceMatrix<=ri;
haptoNbrs(objInd,:) = false; % no haptotaxis towards own nodes
% check which nodes are overlapping
overlapNbrs = distanceMatrix<r_overlap;
overlapNbrs(objInd,:) = false; % allow self-overlap
haptoNbrs(overlapNbrs(:)) = false; % exclude overlapping nodes from haptotaxis

if any(haptoNbrs(:))
    % find unit vectors pointing from neighbours to node
    e_nN = bsxfun(@rdivide,[distanceMatrixFull(haptoNbrs(:)) distanceMatrixFull(find(haptoNbrs(:)) + N*M)],...
        distanceMatrix(haptoNbrs(:))); % bsxfun has similar performace to implicit expansion (below) but is mex-file compatible
    %     e_nN = [distanceMatrixFull(collisionNbrs(:)) distanceMatrixFull(find(collisionNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
    %         ./distanceMatrix(collisionNbrs(:)); % normalise for distance
    Fh = f_hapt.*mean(-e_nN); % - unit vector as we want direction pointing FROM node TO Nbrs
    if any(overlapNbrs(:))
        e_nC = bsxfun(@rdivide,[distanceMatrixFull(overlapNbrs(:)) distanceMatrixFull(find(overlapNbrs(:)) + N*M)],...
            distanceMatrix(overlapNbrs(:))); % directions for overlapping Nbrs
        Fh = Fh + f_hapt.*mean(e_nC); % different sign as force is now repulsive
    end
else
    Fh = [0, 0];
end
end

