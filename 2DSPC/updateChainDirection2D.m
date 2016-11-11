function arrayOut = updateChainDirection2D(arrayNow,arrayPrev,L,rc,bc)
% updates object directions according to update rules 

% issues/to-do's:
% - does it matter if we have discontinuities in our force laws curves?
% - could also implement collisions as stopping or random direction
% - periodic boundaries are only implemented for L>2*rc (I think)
% - mixed periodic boundary conditions can be quite slow
% - extend to enforce segment length
% - move collision detection into a separate function?

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

N = size(arrayPrev,1);
M = size(arrayPrev,2);

% find distances between all pairs of objects
distanceMatrixXY = computeChainDistancesWithBCs(arrayPrev(:,:,[x y]),L,bc);

for objCtr = 1:N
    % calculate force contributions
    
    % motile
    Fm = NaN(2,M);
    % calculate average direction of all nodes in object
    Fm(:,1) = [cos(arrayPrev(objCtr,1,phi)); sin(arrayPrev(objCtr,1,phi))]; % self-alignment (persistence)
    for nodeCtr = 2:M % WARNING this doesn't currently work for periodic boundaries, as the direction can point to the other side of the domain
        Fm(:,nodeCtr) = arrayPrev(objCtr,nodeCtr - 1,[x y]) - arrayPrev(objCtr,nodeCtr,[x y]);% move towards previous node's position
    end
    
    % core repulsion (volume exclusion)
    Fc = NaN(2,M);
    for nodeCtr = 1:M
        Fc(:,nodeCtr) = exclusionForce(squeeze(distanceMatrixXY(:,nodeCtr,:,nodeCtr,:)),objCtr, rc);
    end
    
    % sum forces
    F = Fm + 1e100*Fc;
    % 1e100 to represent Inf vector, but still have direction, should work as long as Fc~O(L)<1.3408e+54
    
    % update directions
    arrayNow(objCtr,:,phi) = atan2(F(y,:),F(x,:));
end

arrayOut = arrayNow;