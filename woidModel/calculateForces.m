function forceArray = calculateForces(arrayPrev,rc,distanceMatrixXY,distanceMatrix,theta,reversals,segmentLength,v_target)
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

distanceMatrixXY = permute(distanceMatrixXY,[3 4 5 1 2]); % this will make indexing later on faster without need for squeeze()
% format of distanceMatrixXY is now N by M by [x y] by N by M

% % calculate average direction of all nodes in object
% phi_CoM = mean(arrayPrev(:,:,phi),2);

forceArray = zeros(N,M,2); % preallocate forces

for objCtr = 1:N
    % check if worm is currently reversing
    if ~reversals(objCtr,2)
        headInd = 1;
        bodyInd = 2:M;
    else
        headInd = M;
        bodyInd = (M-1):-1:1;
    end
    movState = 1 - 2*reversals(objCtr,2); % =-1 if worm is reversing, 1 if not
    % calculate force contributions
    F = NaN(M,2);
    
    % motile
    Fm = NaN(M,2);
    angle = wrapToPi(arrayPrev(objCtr,headInd,phi) ... % previous direction
    + diff(theta(objCtr,headInd,:)) ...% change in internal oscillator
    + pi*diff(reversals(objCtr,:))); % 180 degree turn when reversal starts or ends
    Fm(headInd,:) = [cos(angle), sin(angle)]; 
    for nodeCtr = bodyInd % WARNING this doesn't currently work for periodic boundaries, as the direction can point to the other side of the domain
        Fm(nodeCtr,:) = arrayPrev(objCtr,nodeCtr - 1*movState,[x y]) ...
            - arrayPrev(objCtr,nodeCtr,[x y]);% move towards previous node's position    
    end
    
    angles = atan2(Fm(bodyInd,y),Fm(bodyInd,x)) - diff(theta(objCtr,bodyInd,:),1,3)'; % undulations incl phase shift along worm
    Fm(bodyInd,:) = [cos(angles), sin(angles)];
    % fix magnitue of motile force to give target velocity
    Fm = v_target(objCtr).*Fm;
    
%     % core repulsion (volume exclusion)
%     Fc = NaN(2,M);
%     for nodeCtr = 1:M
%         Fc(:,nodeCtr) = exclusionForce(distanceMatrixXY(:,:,:,objCtr,nodeCtr),...
%             distanceMatrix(:,:,objCtr,nodeCtr),objCtr,nodeCtr,2*rc); % factor of two so that rc is node radius
%     end
    
    % length constraint
    dl = squeeze(arrayPrev(objCtr,2:M,[x y]) - arrayPrev(objCtr,1:M-1,[x y])) ... % direction to next node
        .*repmat(diag(squeeze(distanceMatrix(objCtr,2:M,objCtr,1:M-1))) - segmentLength,1,2); % deviation from segmentLength
    Fl = 4.*([dl; 0 0] - [0 0; dl]); % add forces to next and previous nodes shifted
    
%     % sum motile and exclusion forces with equal magnitude
%     for nodeCtr = 1:M
%         if any(Fc(:,nodeCtr))
%             F(:,nodeCtr) = Fm(:,nodeCtr)./sqrt(sum(Fm(:,nodeCtr).^2)) ...
%                 + Fc(:,nodeCtr)./sqrt(sum(Fc(:,nodeCtr).^2));
%         else
%             F(:,nodeCtr) = Fm(:,nodeCtr);
%         end
%     end
    F = Fm + Fl;
    
    forceArray(objCtr,:,:) = F;
end