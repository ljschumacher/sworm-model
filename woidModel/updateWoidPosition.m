function arrayOut = updateWoidPosition(arrayNow,arrayPrev,v0,bc,L,segmentLength)
% update positions based on current directions

% issues/to-do's:
%   - assert length constraints enforcement
%   - length can be violated again by checking boundary conditions (though
%   should be met considering arclength of high M limit)

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

M = size(arrayPrev,2); % number of nodes

% update position
arrayNow(:,:,x) = arrayPrev(:,:,x) + ...
    v0*cos(arrayNow(:,:,phi));
arrayNow(:,:,y) = arrayPrev(:,:,y) + ...
    v0*sin(arrayNow(:,:,phi));

% enforce length constraints within Woid
% reset positions to one segmentLength away from previous node, along line
% of connection - % WARNING this doesn't currently work for periodic boundaries, as the direction can point to the other side of the domain
for nodeCtr = 2:M
    segmentVec = arrayNow(:,nodeCtr,[x y]) - arrayNow(:,nodeCtr - 1,[x y]); % vec from prev to current node's current pos
    arrayNow(:,nodeCtr,[x y]) = arrayNow(:,nodeCtr - 1,[x y]) ...%prev node's current pos
        + segmentLength*(segmentVec)./...
        repmat(sqrt(sum(segmentVec.^2,3)),1,1,2); % normalise for length of vec connecting nodes
end

arrayNow = checkWoidBoundaryConditions(arrayNow,bc,L);

arrayOut = arrayNow;