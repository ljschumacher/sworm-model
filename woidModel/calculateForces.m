function forceArray = calculateForces(arrayPrev,rc,distanceMatrixXY,distanceMatrix,...
    theta,reversals,segmentLength,v_target,k_l)
% updates object directions according to update rules

% issues/to-do's:
% - mixed periodic boundary conditions can be quite slow
% - calculate forces without loops?

% short-hand for indexing coordinates
x =     1; 
y =     2;
phi =   3;

N = size(arrayPrev,1);
M = size(arrayPrev,2);

distanceMatrixXY = permute(distanceMatrixXY,[3 4 5 1 2]); % this will make indexing later on faster without need for squeeze()
% format of distanceMatrixXY is now N by M by [x y] by N by M

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
    % length constraint
    dl = squeeze(arrayPrev(objCtr,2:M,[x y]) - arrayPrev(objCtr,1:M-1,[x y])) ... % direction to next node
        .*repmat(diag(squeeze(distanceMatrix(objCtr,2:M,objCtr,1:M-1))) - segmentLength,1,2); % deviation from segmentLength
    Fl = k_l.*([dl; 0 0] - [0 0; dl]); % add forces to next and previous nodes shifted
    % sum force contributions
    forceArray(objCtr,:,:) = Fm + Fl;
end
% resolve contact forces
Fc = NaN(N,M,2);
for objCtr = 1:N
    for nodeCtr = 1:M
        Fc(objCtr,nodeCtr,:) = resolveContacts(forceArray,distanceMatrixXY(:,:,:,objCtr,nodeCtr),...
            distanceMatrix(:,:,objCtr,nodeCtr),objCtr,nodeCtr,2*rc); % factor of two so that rc is node radius
    end
end
forceArray = forceArray + Fc;