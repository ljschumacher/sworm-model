function [xyOut, headingsOut] = applyForces(arrayPrev,forceArray,dT,headings,bc,L)
% update positions based on current directions, respecting boundary
% conditions

% issues/to-do's:
%   - length can be violated again by boundary conditions (though
%   should be met considering arclength of high M limit)

% short-hand for indexing coordinates
x =     1;
y =     2;

% N = size(arrayPrev,1); % number of woids
% M = size(arrayPrev,2); % number of nodes

arrayNow = arrayPrev;

% forceAngles = atan2(forceArray(:,:,y),forceArray(:,:,x));
v = sqrt(sum(forceArray.^2,3));

% update position - forward Euler
arrayNow(:,:,x) = arrayPrev(:,:,x) + forceArray(:,:,x)*dT;
arrayNow(:,:,y) = arrayPrev(:,:,y) + forceArray(:,:,y)*dT;

% assert no large displacements - max(v)*dT = v0*dT0
assert(~any(abs(arrayPrev(:) - arrayNow(:))>4*max(v(:))*dT),...
    'Uh-oh, something has gone wrong... (large displacements)')

% correct heading (e.g. if movement has been constrained)
headings = correctHeading(forceArray,headings,bc,L);

[xyOut, headingsOut] = checkWoidBoundaryConditions(arrayNow,headings,bc,L);
