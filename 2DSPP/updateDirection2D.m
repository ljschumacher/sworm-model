function arrayOut = updateDirection2D(arrayNow,arrayPrev,L,rc,bc)
% updates object directions according to update rules 

% issues/to-do's:
% - does it matter if we have discontinuities in our force laws curves?
% - could also implement collisions as stopping or random direction
% - periodic boundaries are only implemented for L>2*rc (I think)
% - mixed periodic boundary conditions can be quite slow
% - how to add in a motile force as well as alignment and exclusion?
% - extend to objects made of multiple nodes
% - move collision detection into a separate function?

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

N = size(arrayPrev,1);
ndim = 2;

% find distances between all pairs of objects
distanceMatrixXY = computeDistancesWithBCs(arrayPrev(:,[x y]),L,bc);
distanceMatrix = sqrt(sum(distanceMatrixXY.^2,3)); % reduce to scalar


for objCtr = 1:N
    % calculate force contributions
    
    % alignment
    % calculate average direction of all nodes in object
    Fa = [cos(arrayPrev(objCtr,phi)); sin(arrayPrev(objCtr,phi))]; % self-alignment (persistence)
    
    % core repulsion (volume exclusion)
    collisionNbrs = distanceMatrix(:,objCtr)<=rc;
    collisionNbrs(objCtr) = false; % no self-repulsion
    Nc = nnz(collisionNbrs);
    Fc = sum(... % sum over all collision neighbours
        reshape(distanceMatrixXY(objCtr,collisionNbrs,:),Nc,ndim)... %direction FROM neighbours TO object (use reshape rather than squeeze for case when Nc = 1)
        ./distanceMatrix(collisionNbrs,objCtr*ones(1,ndim))... % normalise for distance
        ,1)'; % made into column vector
    
    % sum forces
    F = Fa + 1e100*Fc;
    % 1e100 to represent Inf vector, but still have direction, should work as long as Fc~O(L)<1.3408e+54
    
    % update directions
    arrayNow(objCtr,phi) = atan2(F(y),F(x));
end

arrayOut = arrayNow;