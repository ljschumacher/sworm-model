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

forceAngles = atan2(forceArray(:,:,y),forceArray(:,:,x));
v = sqrt(sum(forceArray.^2,3));
assert(~any(isinf(v(:))),'In this code we strive to respect the laws of physics...')

% update position
arrayNow(:,:,x) = arrayPrev(:,:,x) + ...
    v.*cos(forceAngles);
arrayNow(:,:,y) = arrayPrev(:,:,y) + ...
    v.*sin(forceAngles);

% update direction
arrayNow(:,:,phi) = forceAngles; % not sure if this is the right use of phi, which technically points along the connection btw nodes

arrayOut = checkWoidBoundaryConditions(arrayNow,bc,L);
