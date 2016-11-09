function [outArray] = initialiseObjects2D(inArray,L)
% initialises object positions and directions
% uniformly randomly distributed

% issues/to-do's:
% - initial positions do not currently respect volume exclusion
% - could implement non-random initial positions, e.g. regularly spaced

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

N = size(inArray,1);
if size(L,1)>size(L,2)
    L = L';
end

for objCtr = 1:N
    % Position, should work for both scalar and vector L
    inArray(objCtr,[x y],1) = L.*rand(1,2);
    
    % Direction
    inArray(objCtr,phi,1) = pi*(2*rand - 1);   % phi between -pi and pi
end

outArray = inArray;

