function [outArray] = initialiseWoids(inArray,L,segmentLength,theta,rc,bc)
% initialises object positions and directions
% uniformly randomly distributed

% issues/to-do's:
% - needlessy re-computing ditances between all successfully initialised
% agents when checking for overlaps
% - not tested for periodic boundaries


% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

N = size(inArray,1);
M = size(inArray,2);
if size(L,1)>size(L,2)
    L = L';
end

L0 = (M-1)*segmentLength;
disp('Initialising woid positions...')
objCtr = 1;
while objCtr <= N
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
    % check for any overlaps, and if so discard this woid's position
    distanceMatrix = sqrt(sum(computeWoidDistancesWithBCs(inArray(1:objCtr,:,[x y],1),L,bc).^2,5));% scalar distance from current object to any other
    collisionNbrs = squeeze(any(any(2*rc>=distanceMatrix(objCtr,:,:,:),4),2));
    collisionNbrs(objCtr) = false; % don't volume-exclude self
    if ~any(collisionNbrs) 
        objCtr = objCtr + 1;
    end
end

outArray = inArray;
