function arrayOut = updateWoidDirection(arrayNow,arrayPrev,L,rc,bc,theta)
% updates object directions according to update rules

% issues/to-do's:
% - does it matter if we have discontinuities in our force laws curves?
% - mixed periodic boundary conditions can be quite slow
% - improve loop structure and motile vs exclusion force

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

N = size(arrayPrev,1);
M = size(arrayPrev,2);

% find distances between all pairs of objects
distanceMatrixXY = computeWoidDistancesWithBCs(arrayPrev(:,:,[x y]),L,bc);
distanceMatrixXY = permute(distanceMatrixXY,[3 4 5 1 2]); % this will make indexing later on faster without need for squeeze()

% % calculate average direction of all nodes in object
% phi_CoM = mean(arrayPrev(:,:,phi),2);

for objCtr = 1:N
    % calculate force contributions
    F = NaN(2,M);
    
    % motile
    Fm = NaN(2,M);
    angle = arrayPrev(objCtr,1,phi) + diff(theta(objCtr,1,:));
    Fm(:,1) = [cos(angle); sin(angle)]; 
    for nodeCtr = 2:M % WARNING this doesn't currently work for periodic boundaries, as the direction can point to the other side of the domain
        Fm(:,nodeCtr) = arrayPrev(objCtr,nodeCtr - 1,[x y]) ...
            - arrayPrev(objCtr,nodeCtr,[x y]);% move towards previous node's position    
    end
    angles = arrayPrev(objCtr,2:end,phi) + diff(theta(objCtr,2:end,:),1,3); % undulations incl phase shift along worm
    Fm(:,2:end) = [cos(angles); sin(angles)];
    
    % core repulsion (volume exclusion)
    Fc = NaN(2,M);
    for nodeCtr = 1:M
        Fc(:,nodeCtr) = exclusionForce(distanceMatrixXY(:,:,:,objCtr,nodeCtr),...
            objCtr,nodeCtr,2*rc); % factor of two so that rc is node radius
    end
    
    % sum motile and exclusion forces with equal magnitude
    for nodeCtr = 1:M
        if any(Fc(:,nodeCtr))
            F(:,nodeCtr) = Fm(:,nodeCtr)./sqrt(sum(Fm(:,nodeCtr).^2)) ...
                + Fc(:,nodeCtr)./sqrt(sum(Fc(:,nodeCtr).^2));
        else
            F(:,nodeCtr) = Fm(:,nodeCtr);
        end
    end
    
    % update directions
    arrayNow(objCtr,:,phi) = atan2(F(y,:),F(x,:));
end

arrayOut = arrayNow;