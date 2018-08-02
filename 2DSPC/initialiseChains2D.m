function [outArray] = initialiseChains2D(inArray,L,segmentLength,deltaTheta)
% initialises object positions and directions
% uniformly randomly distributed

% issues/to-do's:
% - initial positions do not currently respect volume exclusion or chain
% overlap


% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

N = size(inArray,1);
M = size(inArray,2);
if size(L,1)>size(L,2)
    L = L';
end

L0 = M*segmentLength;

for objCtr = 1:N
    % initialise head node
    % Position, should work for both scalar and vector L
    inArray(objCtr,1,[x y],1) = L0 + (L - 2*L0).*rand(1,2); % initialise positions at least one chain length away from edge
    % Direction
    inArray(objCtr,1,phi,1) = pi*(2*rand - 1);   % phi between -pi and pi
    for nodeCtr = 2:M % initialise chain positions, node by node
        % initialise the next node in the right direction at segmentLength away
        inArray(objCtr,nodeCtr,phi,1) = wrapToPi(inArray(objCtr,nodeCtr - 1,phi,1) ...% previous node's orientation
            + deltaTheta*(2*rand - 1));% choose from within ± some angle of previous node's orientation
        inArray(objCtr,nodeCtr,[x y],1) = squeeze(inArray(objCtr,nodeCtr - 1,[x y],1))... % previous node's position
            - segmentLength*[cos(inArray(objCtr,nodeCtr,phi,1)); sin(inArray(objCtr,nodeCtr,phi,1))];
    end
end

outArray = inArray;
