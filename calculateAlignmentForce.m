function [Fa] = calculateAlignmentForce(headings,distanceMatrix,objInd,headInd,ri,f_align,selfAlign)
% calculates the alignment contribution to the motile force, i.e. the
% mean direction of Nbrs within an interaction radius
N = size(distanceMatrix,1);
M = size(distanceMatrix,2);
% check which nodes are within interaction radius
alignNbrs = distanceMatrix<=ri;
alignNbrs(objInd,:) = false; % no alignment towards own (non-head) nodes
alignNbrs(objInd,headInd) = selfAlign; % alignment towards own direction

if any(alignNbrs(:))
%     % find unit vectors pointing from neighbours to node
%     e_nN = bsxfun(@rdivide,[distanceMatrixFull(alignNbrs(:)) distanceMatrixFull(find(alignNbrs(:)) + N*M)],...
%         distanceMatrix(alignNbrs(:))); % bsxfun has similar performace to implicit expansion (below) but is mex-file compatible
%     %     e_nN = [distanceMatrixFull(collisionNbrs(:)) distanceMatrixFull(find(collisionNbrs(:)) + N*M)]... %direction FROM neighbours TO object [x, y]
%     %         ./distanceMatrix(collisionNbrs(:)); % normalise for distance
    meanHeading = mean(headings(alignNbrs(:)));
    Fa = f_align.*[cos(meanHeading), sin(meanHeading)];
else
    Fa = [0, 0];
end
end

