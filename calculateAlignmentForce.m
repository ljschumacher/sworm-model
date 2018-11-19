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
    meanHeading = mean(headings(alignNbrs(:)));
    Fa = f_align.*[cos(meanHeading), sin(meanHeading)];
else
    Fa = [0, 0];
end
end

