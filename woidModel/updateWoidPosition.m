function arrayOut = updateWoidPosition(arrayNow,arrayPrev,v,bc,L,segmentLength,reversals)
% -- DEPRECATED --

% update positions based on current directions

% issues/to-do's:
%   - assert length constraints enforcement
%   - length can be violated again by boundary conditions (though
%   should be met considering arclength of high M limit)
%   - could move length constraint as force into updateWoidDirections or
%   separate function

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

N = size(arrayPrev,1); % number of woids
M = size(arrayPrev,2); % number of nodes

% update position
if length(v)==1
    arrayNow(:,:,x) = arrayPrev(:,:,x) + ...
        v*cos(arrayNow(:,:,phi));
    arrayNow(:,:,y) = arrayPrev(:,:,y) + ...
        v*sin(arrayNow(:,:,phi));
elseif length(v)==N
    arrayNow(:,:,x) = arrayPrev(:,:,x) + ...
        repmat(v,1,M).*cos(arrayNow(:,:,phi));
    arrayNow(:,:,y) = arrayPrev(:,:,y) + ...
        repmat(v,1,M).*sin(arrayNow(:,:,phi));
else
    error('Number of elements in velocity vectors must be 1 or N')
end

% % enforce length constraints within Woid
% % reset positions to one segmentLength away from previous node, along line
% % of connection - % WARNING this doesn't currently work for periodic boundaries, as the direction can point to the other side of the domain
% fwdInd = find(~reversals); % all the woids that are moving forward
% if any(fwdInd)
%     for nodeCtr = 2:M
%         segmentVec = arrayNow(fwdInd,nodeCtr,[x y]) - arrayNow(fwdInd,nodeCtr - 1,[x y]); % vec from prev to current node's current pos
%         arrayNow(fwdInd,nodeCtr,[x y]) = arrayNow(fwdInd,nodeCtr - 1,[x y]) ...%prev node's current pos
%             + segmentLength*(segmentVec)./...
%             repmat(sqrt(sum(segmentVec.^2,3)),1,1,2); % normalise for length of vec connecting nodes
%         % %     % update direction to how the node actually moved
%         % %     dx = arrayNow(:,nodeCtr,[x y]) - arrayPrev(:,nodeCtr,[x y]);
%         % %     arrayNow(:,nodeCtr,phi) = atan2(dx(:,y),dx(:,x));
%     end
% end
% bwdInd = find(reversals); % all the woids that are moving forward
% if any(bwdInd)
%     for nodeCtr = (M-1):-1:1
%         segmentVec = arrayNow(bwdInd,nodeCtr,[x y]) - arrayNow(bwdInd,nodeCtr + 1,[x y]); % vec from prev to current node's current pos
%         arrayNow(bwdInd,nodeCtr,[x y]) = arrayNow(bwdInd,nodeCtr + 1,[x y]) ...%prev node's current pos
%             + segmentLength*(segmentVec)./...
%             repmat(sqrt(sum(segmentVec.^2,3)),1,1,2); % normalise for length of vec connecting nodes
%     end
% end
arrayNow = checkWoidBoundaryConditions(arrayNow,bc,L);

arrayOut = arrayNow;