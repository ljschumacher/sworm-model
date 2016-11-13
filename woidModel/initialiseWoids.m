function [outArray] = initialiseWoids(inArray,L,segmentLength,theta)
% initialises object positions and directions
% uniformly randomly distributed

% issues/to-do's:
% - initial positions do not currently respect volume exclusion 


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
    inArray(objCtr,1,[x y],1) = L0 + (L - 2*L0).*rand(1,2); % initialise positions at least one woid length away from edge
    % Direction
    inArray(objCtr,:,phi,1) = wrapToPi(pi*(2*rand - 1) - theta(objCtr,:));   % random orientation between -pi and pi for each object plus undulations
    for nodeCtr = 2:M % initialise woid positions, node by node
        % initialise the next node in the right direction at segmentLength away
        inArray(objCtr,nodeCtr,[x y],1) = squeeze(inArray(objCtr,nodeCtr - 1,[x y],1))... % previous node's position
            - segmentLength*[cos(inArray(objCtr,nodeCtr,phi,1)); sin(inArray(objCtr,nodeCtr,phi,1))];
    end
end

outArray = inArray;
