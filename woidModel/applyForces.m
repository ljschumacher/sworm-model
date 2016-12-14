function arrayOut = applyForces(arrayPrev,forceArray,bc,L)
% update positions based on current directions, respecting boundary
% conditions

% issues/to-do's:
%   - length can be violated again by boundary conditions (though
%   should be met considering arclength of high M limit)

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

% N = size(arrayPrev,1); % number of woids
% M = size(arrayPrev,2); % number of nodes

arrayNow = arrayPrev;

% update direction
arrayNow(:,:,phi) = atan2(forceArray(:,:,y),forceArray(:,:,x)); % not sure if this is the right use of phi, which technically points along the connection btw nodes
v = sqrt(sum(forceArray.^2,3));

% update position
arrayNow(:,:,x) = arrayPrev(:,:,x) + ...
    v.*cos(arrayNow(:,:,phi));
arrayNow(:,:,y) = arrayPrev(:,:,y) + ...
    v.*sin(arrayNow(:,:,phi));

arrayOut = checkWoidBoundaryConditions(arrayNow,bc,L);
