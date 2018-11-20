function arrayOut = updateWoidPosition(arrayNow,arrayPrev,v,bc,L)
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

arrayNow = checkWoidBoundaryConditions(arrayNow,bc,L);

arrayOut = arrayNow;