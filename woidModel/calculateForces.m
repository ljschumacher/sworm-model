function forceArray = calculateForces(posPrev,posPrevPrev,rc,distanceMatrixXY,distanceMatrix,...
    theta,reversals,segmentLength,v_target,k_l,k_theta, r_LJcutoff, eps_LJ)
% updates object directions according to update rules

% issues/to-do's:
% - mixed periodic boundary conditions can be quite slow
% - calculate forces without loops?
% - once we have the bending constraints, do we still need the undulations
% at every node??
% - refactor individual forces into their own functions
% - could improve tangent estimate for motile force by estimating direction
% towards point on curve
% - should body angle be calculated from improved tangent estimate (like
% motile force?)
% - should the target angles for bending actually be
% theta(objCtr,bodyInd,end)?

% short-hand for indexing coordinates
x =     1;
y =     2;

N = size(posPrev,1);
M = size(posPrev,2);

distanceMatrixXY = permute(distanceMatrixXY,[3 4 5 1 2]); % this will make indexing later on faster without need for squeeze()
% format of distanceMatrixXY is now N by M by [x y] by N by M

forceArray = zeros(N,M,2); % preallocate forces

for objCtr = 1:N
    % check if worm is currently reversing
    if ~reversals(objCtr)
        headInd = 1;
        bodyInd = 2:M;
    else
        headInd = M;
        bodyInd = (M-1):-1:1;
    end
    movState = 1 - 2*reversals(objCtr); % =-1 if worm is reversing, 1 if not
    % calculate force contributions
    % motile force
    Fm = NaN(M,2);
    ds = NaN(M,2); %change in node positon along worm
    % estimate tangent vectors
    % head tangent: direction from next node's position
    ds(headInd,:) = posPrev(objCtr,headInd,[x y]) ...
        - posPrev(objCtr,headInd + 1*movState,[x y]);
    % calculate direction towards previous node's position
    ds(bodyInd,:) = posPrev(objCtr,bodyInd - 1*movState,[x y]) ...
        - posPrev(objCtr,bodyInd,[x y]);
    % before estimating tangents, calculate polar unit vectors
    phi = atan2(ds(bodyInd,y),ds(bodyInd,x));
    e_phi = [-sin(phi) cos(phi)]; % unit vector in direction of phi, size M-1 by 2
    % body tangents: estimate tangent direction by average of directions towards previous node and from next node
    ds(bodyInd(1:end-1),:) = (ds(bodyInd(1:end-1),:) + ds(bodyInd(2:end),:))/2;
    % (tail tangent: towards previous node)
    
    % head motile force
    headAngle = wrapToPi(theta(objCtr,headInd,1) ... % previous heading
        + diff(theta(objCtr,headInd,:),1,3));% change in internal oscillator
    
    Fm(headInd,:) = [cos(headAngle), sin(headAngle)];
    % body motile force
    Fm(bodyInd,:) = ds(bodyInd,:);
    % fix magnitue of motile force to give target velocity
    Fm = v_target(objCtr).*Fm./sqrt(sum(Fm.^2,2));
    
    % length constraint
    dl = squeeze(posPrev(objCtr,2:M,[x y]) - posPrev(objCtr,1:M-1,[x y])) ... % direction to next node
        .*(diag(squeeze(distanceMatrix(objCtr,2:M,objCtr,1:M-1))) - segmentLength) ...% deviation from segmentLength
        ./sqrt(sum(squeeze(posPrev(objCtr,2:M,[x y]) - posPrev(objCtr,1:M-1,[x y])).^2,2)); % normalised for segment length
    Fl = k_l.*([dl; 0 0] - [0 0; dl]); % add forces to next and previous nodes shifted
    
    % bending constraints - rotational springs with changing 'rest length' due to active undulations
    bodyAngles = unwrap(atan2(ds([headInd, bodyInd],y),ds([headInd, bodyInd],x)));
    % previous bodyAngles
    dsPrev(headInd,:) = posPrevPrev(objCtr,headInd,[x y]) ...
        - posPrevPrev(objCtr,headInd + 1*movState,[x y]);
    dsPrev(bodyInd,:) = posPrevPrev(objCtr,bodyInd - 1*movState,[x y]) ...
        - posPrevPrev(objCtr,bodyInd,[x y]);
    dsPrev(bodyInd(1:end-1),:) = (dsPrev(bodyInd(1:end-1),:) + dsPrev(bodyInd(2:end),:))/2;
    bodyAnglesPrev = unwrap(atan2(dsPrev([headInd, bodyInd],y),dsPrev([headInd, bodyInd],x)));
    torques = k_theta.*wrapToPi(gradient(bodyAngles) - gradient(bodyAnglesPrev));
%     torques = k_theta.*wrapToPi(gradient(bodyAngles) - gradient(unwrap(theta(objCtr,[headInd, bodyInd],1)')));
    torques = torques(2:end-1); % no rotational springs at head and tail node
%         bodyAngles = atan2(ds(bodyInd,y),ds(bodyInd,x));
%         torques = k_theta.*wrapToPi(diff(unwrap(bodyAngles))); % this straightens worms
    
    l = sqrt(sum(ds(bodyInd,:).^2,2)); % length between node and prev node, length M-1
    momentsfwd = torques.*l(1:end-1).*e_phi(1:end-1,:);
    momentsbwd = torques.*l(2:end).*e_phi(2:end,:);
    F_theta = NaN(M,2); % pre-allocate to index nodes in order depending on movement state
    F_theta([headInd, bodyInd],:) = [momentsfwd; 0 0; 0 0] ... % rotational force from node n+1 onto n
        + [0 0; 0 0; momentsbwd] ...% rotational force from node n-1 onto n
        + [0 0; -(momentsfwd + momentsbwd); 0 0];% reactive force on node n (balancing forces exerted onto nodes n+1 and n -1
    % sum force contributions
    forceArray(objCtr,:,:) = Fm + Fl + F_theta;
%     % uncomment for debugging...
%     plot(squeeze(posPrev(objCtr,:,x)),squeeze(posPrev(objCtr,:,y)),'.-'), axis equal, hold on
%     quiver(squeeze(posPrev(objCtr,:,x))',squeeze(posPrev(objCtr,:,y))',F_theta(:,1),F_theta(:,2),1)
%     1;
end
% resolve contact forces
Fc = NaN(N,M,2);
for objCtr = 1:N
    for nodeCtr = 1:M
        Fc(objCtr,nodeCtr,:) = resolveContacts(forceArray,distanceMatrixXY(:,:,:,objCtr,nodeCtr),...
            distanceMatrix(:,:,objCtr,nodeCtr),objCtr,nodeCtr,2*rc,r_LJcutoff,eps_LJ); % factor of two so that rc is node radius
    end
end
forceArray = forceArray + Fc;