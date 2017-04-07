function [xyOut, thetaOut] = applyForces(arrayPrev,forceArray,theta,bc,L)
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

forceAngles = atan2(forceArray(:,:,y),forceArray(:,:,x));
v = sqrt(sum(forceArray.^2,3));

% update position
arrayNow(:,:,x) = arrayPrev(:,:,x) + ...
    v.*cos(forceAngles);
arrayNow(:,:,y) = arrayPrev(:,:,y) + ...
    v.*sin(forceAngles);

% correct heading (e.g. if movement has been constrained)
theta = correctHeading(forceArray,theta,bc,L);
    
[xyOut, thetaOut] = checkWoidBoundaryConditions(arrayNow,theta,bc,L);
