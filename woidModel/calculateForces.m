function forceArray = calculateForces(distanceMatrixXY,distanceMatrix,rc,...
    theta,reversals,segmentLength,v_target,k_l,k_theta,phaseOffset,sigma_LJ,r_LJcutoff, eps_LJ)
% updates object directions according to update rules

% issues/to-do's:
% - mixed periodic boundary conditions can be quite slow
% - calculate forces without loops?
% - refactor individual forces into their own functions
% - M>1 forces need fixing for periodic boundaries - vectorial directions
% aren't corrected (could be fixed by calculating from distanceMatrixXY?
% Which may also enable vectorization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%

% short-hand for indexing coordinates
x =     1;
y =     2;

N = size(distanceMatrixXY,1);
M = size(distanceMatrixXY,2);

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
    if M>1
        % head tangent: direction from next node's position
        ds(1,:) = distanceMatrixXY(objCtr,2,[x y],objCtr,1);
        % calculate direction towards previous node's position
        ds(2:M,:) = distanceMatrixXY(objCtr,2:M,[x y],objCtr,1:M-1);
        % before estimating tangents, calculate polar unit vectors
        phi = atan2(ds(2:M,y),ds(2:M,x));
        e_phi = [-sin(phi) cos(phi)]; % unit vector in direction of phi, size M-1 by 2
        l = sqrt(sum(ds(2:M,:).^2,2)); % length between node and prev node, length M-1
        % body tangents: estimate tangent direction by average of directions towards previous node and from next node
        ds(2:M-1,:) = (ds(2:M-1,:) + ds(3:M,:))/2;
        % (tail tangent: towards previous node)
    end
    
    % head motile force
    headAngle = theta(objCtr,headInd);
    
    Fm(headInd,:) = [cos(headAngle), sin(headAngle)];
    % body motile force
    Fm(bodyInd,:) = movState*ds(bodyInd,:);
    % fix magnitue of motile force to give target velocity
    Fm = v_target(objCtr).*Fm./sqrt(sum(Fm.^2,2));
    
    % length constraint
    if k_l>0&&M>1
        dl = distanceMatrixXY(objCtr,1:M-1,[x y],objCtr,2:M)./l' ... % direction to next node, normalised for segment length
            .*(l' - segmentLength);% deviation from segmentLength
        nl = 1./(1 - sum(dl.^2,3)/segmentLength.^2); % non-linear part of spring
        Fl = k_l.*squeeze(cat(2,dl.*nl,zeros(1,1,2)) - cat(2,zeros(1,1,2),dl.*nl)); % add forces to next and previous nodes shifted
    else
        Fl = zeros(size(Fm));
    end
    
    % bending constraints - rotational springs with changing rest angle due to active undulations
    if k_theta(objCtr)>0&&M>2
        bodyAngles = unwrap(atan2(ds(:,y),ds(:,x)));
        dK = (gradient(bodyAngles) - gradient(sin(phaseOffset(objCtr,:))'));
        %         %     % uncomment for debugging...
        %     plot(gradient(bodyAngles)), hold on, plot(gradient(sin(phaseOffset(objCtr,:))'))
        torques = k_theta(objCtr).*dK;%./(1 - dK.^2/pi^2);
        torques = torques(2:end-1); % no rotational springs at head and tail node
        
        momentsfwd = torques.*l(1:end-1).*e_phi(1:end-1,:);
        momentsbwd = torques.*l(2:end).*e_phi(2:end,:);
        %     F_theta = NaN(M,2); % pre-allocate to index nodes in order depending on movement state
        F_theta = [momentsfwd; 0 0; 0 0] ... % rotational force from node n+1 onto n
            + [0 0; 0 0; momentsbwd] ...% rotational force from node n-1 onto n
            + [0 0; -(momentsfwd + momentsbwd); 0 0];% reactive force on node n (balancing forces exerted onto nodes n+1 and n -1
    else
        F_theta = zeros(size(Fm));
    end
    % sum force contributions
    forceArray(objCtr,:,:) = Fm + Fl + F_theta;
    %     % uncomment for debugging...
    %         plot(squeeze(posPrev(objCtr,:,x)),squeeze(posPrev(objCtr,:,y)),'.-'), axis equal, hold on
    %         quiver(squeeze(posPrev(objCtr,:,x))',squeeze(posPrev(objCtr,:,y))',Fm(:,1),Fm(:,2),0)
    %         1;
end
% resolve contact forces
Fc = NaN(N,M,2);
for objCtr = 1:N
    for nodeCtr = 1:M
        Fc(objCtr,nodeCtr,:) = resolveContacts(forceArray,distanceMatrixXY(:,:,:,objCtr,nodeCtr),...
            distanceMatrix(:,:,objCtr,nodeCtr),objCtr,nodeCtr,2*rc,sigma_LJ,r_LJcutoff,eps_LJ); % factor of two so that rc is node radius
    end
end
forceArray = forceArray + Fc;