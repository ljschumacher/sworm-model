function forceArray = calculateForces(arrayPrev,rc,distanceMatrixXY,distanceMatrix,...
    theta,reversals,segmentLength,v_target,k_l,k_theta)
% updates object directions according to update rules

% issues/to-do's:
% - mixed periodic boundary conditions can be quite slow
% - calculate forces without loops?
% - once we have the bending constraints, do we still need the undulations
% at every node??
% - refactor individual forces into their own functions

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
    % motile force
    Fm = NaN(M,2);
    ds = NaN(M,2); %change in node positon
    % head motile force
    angle = wrapToPi(arrayPrev(objCtr,headInd,phi) ... % previous direction
        + diff(theta(objCtr,headInd,:)) ...% change in internal oscillator
        + pi*diff(reversals(objCtr,:))); % 180 degree turn when reversal starts or ends
    Fm(headInd,:) = [cos(angle), sin(angle)];
    % body motile force
    ds(bodyInd,:) = arrayPrev(objCtr,bodyInd - 1*movState,[x y]) ...
        - arrayPrev(objCtr,bodyInd,[x y]);% direction towards previous node's position
    bodyAngles = atan2(ds(bodyInd,y),ds(bodyInd,x));
    targetAngles = bodyAngles + diff(theta(objCtr,bodyInd,:),1,3)'; % undulations incl phase shift along worm
    Fm(bodyInd,:) = [cos(targetAngles), sin(targetAngles)];
    % fix magnitue of motile force to give target velocity
    Fm = v_target(objCtr).*Fm;
    % length constraint
    dl = squeeze(arrayPrev(objCtr,2:M,[x y]) - arrayPrev(objCtr,1:M-1,[x y])) ... % direction to next node
        .*repmat((diag(squeeze(distanceMatrix(objCtr,2:M,objCtr,1:M-1))) - segmentLength) ...% deviation from segmentLength
        ./sqrt(sum(squeeze(arrayPrev(objCtr,2:M,[x y]) - arrayPrev(objCtr,1:M-1,[x y])).^2,2)),1,2); % normalised for segment length
    Fl = k_l.*([dl; 0 0] - [0 0; dl]); % add forces to next and previous nodes shifted
    % bending constraints - rotational springs with changing 'rest length' due
    % to active undulations
    torques = k_theta.*(wrapToPi(diff(bodyAngles)) - wrapToPi(diff(targetAngles)));
    %     - diff(theta(objCtr,bodyInd,end),1,2)');% deviation from target change in angle to previous node, length M-2
    e_phi = [-sin(bodyAngles) cos(bodyAngles)]; % unit vector in direction of phi, size M-1 by 2
    l = sqrt(sum(ds(bodyInd,:).^2,2)); % length between node and prev node, length M-1
    momentsfwd = repmat(torques.*l(1:end-1),1,2).*e_phi(1:end-1,:);
    momentsbwd = repmat(torques.*l(2:end),1,2).*e_phi(2:end,:);
    F_theta = NaN(M,2); % pre-allocate to index nodes in order depending on movement state
    F_theta([headInd, bodyInd],:) = [momentsfwd; 0 0; 0 0] ... % rotational force from node n+1 onto n
        + [0 0; 0 0; momentsbwd] ...% rotational force from node n-1 onto n
        + [0 0; -(momentsfwd + momentsbwd); 0 0];% reactive force on node n (balancing forces exerted onto nodes n+1 and n -1
    % sum force contributions
    forceArray(objCtr,:,:) = Fm + Fl + F_theta;
    %     % uncomment for debugging...
    %     plot(squeeze(arrayPrev(objCtr,:,x)),squeeze(arrayPrev(objCtr,:,y)),'.-'), axis equal, hold on
    %     quiver(squeeze(arrayPrev(objCtr,:,x))',squeeze(arrayPrev(objCtr,:,y))',Fm(:,1),Fm(:,2),0)
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