function forceArray = calculateForces(distanceMatrixXY,distanceMatrix,rc,...
    headings,reversals,segmentLength,v_target,k_l,k_theta,theta_0,phaseOffset,...
    sigma_LJ,r_LJcutoff,eps_LJ,LJnodes,LJmode,angleNoise,ri,r_overlap,f_hapt,haptotaxisMode,roamingLogInd,f_align)
% updates object directions according to update rules

% issues/to-do's:
% - mixed periodic boundary conditions can be quite slow
% - calculate forces without loops?
% - refactor individual forces into their own functions
% - if using noise for head angle, better use von mises random number?
% - if changing head orientation based on forces instantaneously, move
% angleNoise generation to updateOscillators function

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
        for nodeCtr = 2:M
            ds(nodeCtr,:) = distanceMatrixXY(objCtr,nodeCtr,[x y],objCtr,nodeCtr-1);
        end
        % before estimating tangent vectors, calculate polar unit vectors and
        % bending angles
        phi = atan2(ds(2:M,y),ds(2:M,x));
        e_phi = [-sin(phi) cos(phi)]; % unit vector in direction of phi, size M-1 by 2
        l = sqrt(sum(ds(2:M,:).^2,2)); % length between node and prev node, length M-1
        bendingAngles = unwrap(atan2(ds(:,y),ds(:,x)));
        % body tangents: estimate tangent direction by average of directions towards previous node and from next node
        ds(2:M-1,:) = (ds(2:M-1,:) + ds(3:M,:))/2;
        % (tail tangent: towards previous node)
    end
    
    % head motile force
    headAngle = headings(objCtr,headInd) + angleNoise*randn();
    
    Fm(headInd,:) = [cos(headAngle), sin(headAngle)];
    % haptotaxis - move to the direction of other worms (should this be added before angular noise?)
    if f_hapt~=0 % haptotaxis could be attractive or repulsive
        if strcmp(haptotaxisMode,'weighted_additive')%&&~roamingLogInd(objCtr) % only update those worms which are not roaming
            Fm(headInd,:) = Fm(headInd,:) ...
                + calculateHaptotaxisWeighted(distanceMatrixXY(:,:,:,objCtr,headInd),...
                distanceMatrix(:,:,objCtr,headInd),objCtr,ri,r_overlap,f_hapt);
        elseif strcmp(haptotaxisMode,'weighted')%&&~roamingLogInd(objCtr) % only update those worms which are not roaming
            Fm(headInd,:) = Fm(headInd,:) ...
                + calculateHaptotaxisWeighted(distanceMatrixXY(:,:,:,objCtr,headInd),...
                distanceMatrix(:,:,objCtr,headInd),objCtr,ri,r_overlap,f_hapt);
        elseif strcmp(haptotaxisMode,'constant')
            Fm(headInd,:) = Fm(headInd,:) ...
                + calculateHaptotaxis(distanceMatrixXY(:,:,:,objCtr,headInd),...
                distanceMatrix(:,:,objCtr,headInd),objCtr,ri,r_overlap,f_hapt);
        end
    end
    
    % Vicsek-type alignment force - only really used for demonstration
    if f_align~=0
        Fm(headInd,:) = Fm(headInd,:) ...
            + calculateAlignmentForce(headings,distanceMatrix(:,:,objCtr,headInd),objCtr,headInd,ri,f_align,true);
    end
    
    % body motile force
    Fm(bodyInd,:) = movState*ds(bodyInd,:);
    % fix magnitue of motile force to give target velocity
    Fm = v_target(objCtr).*Fm./sqrt(sum(Fm.^2,2));
    
    % length constraint
    if k_l>0&&M>1
        dl = NaN(M-1,2);
        for nodeCtr = 1:M-1
            dl(nodeCtr,:) = distanceMatrixXY(objCtr,nodeCtr,[x y],objCtr,nodeCtr+1); % direction to next node
        end
        dl = dl./l.*(l - segmentLength);% normalised for segment length and deviation from segmentLength
        dl_max = segmentLength;
        dl_nl = min(abs(l - segmentLength),0.99*dl_max);% non-linear contribution saturates to avoid force decrease for
        % l - segL > segLmax, which could happen for finite timestep
        nl = 1./abs(1 - (dl_nl/dl_max).^2); % non-linear part of spring
        Fl = k_l.*(cat(1,dl.*nl,zeros(1,2)) - cat(1,zeros(1,2),dl.*nl)); % add forces to next and previous nodes shifted
    else
        Fl = zeros(size(Fm));
    end
    
    % bending constraints - rotational springs with changing rest angle due to active undulations
    if k_theta(objCtr)>0&&M>2
        dK_max = pi;
        if theta_0>0
            dK = (diff(bendingAngles) - diff(sin(phaseOffset(objCtr,:))'));
        else
            dK = diff(bendingAngles);
        end
        dK_nl = min(abs(dK),0.99*dK_max);% non-linear contribution saturates to avoid force decrease for
        % dK > dK_max, which could happen for finite timestep
        torques = k_theta(objCtr).*dK./(1 - dK_nl.^2/dK_max^2);
        torques = torques(2:end); % no rotational springs at head and tail node
        %     % uncomment for debugging...
        %     plot(gradient(bodyAngles)), hold on, plot(gradient(sin(phaseOffset(objCtr,:))'))
        
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
end
% resolve contact forces
if N==40&&M==49 % check if we can use compiled mex function
    Fc = resolveContactsLoop_mex(forceArray,distanceMatrixXY,...
        distanceMatrix,2*rc,sigma_LJ,r_LJcutoff,eps_LJ, LJnodes, LJmode);
elseif N==40&&M==36 % check if we can use compiled mex function
    Fc = resolveContactsLoop_N40_M36_mex(forceArray,distanceMatrixXY,...
        distanceMatrix,2*rc,sigma_LJ,r_LJcutoff,eps_LJ, LJnodes, LJmode);
elseif N==40&&M==18 % check if we can use compiled mex function
    Fc = resolveContactsLoop_N40_M18_mex(forceArray,distanceMatrixXY,...
        distanceMatrix,2*rc,sigma_LJ,r_LJcutoff,eps_LJ, LJnodes, LJmode);
elseif N==200&&M==18 % check if we can use compiled mex function
    Fc = resolveContactsLoop_N200_M18_mex(forceArray,distanceMatrixXY,...
        distanceMatrix,2*rc,sigma_LJ,r_LJcutoff,eps_LJ, LJnodes, LJmode);
elseif N==200&&M==36 % check if we can use compiled mex function
    Fc = resolveContactsLoop_N200_M36_mex(forceArray,distanceMatrixXY,...
        distanceMatrix,2*rc,sigma_LJ,r_LJcutoff,eps_LJ, LJnodes, LJmode);
else
    Fc = NaN(N,M,2);
    for objCtr = 1:N
        for nodeCtr = 1:M
            Fc(objCtr,nodeCtr,:) = resolveContacts(forceArray,distanceMatrixXY(:,:,:,objCtr,nodeCtr),...
                distanceMatrix(:,:,objCtr,nodeCtr),objCtr,nodeCtr,2*rc,...
                sigma_LJ,r_LJcutoff,eps_LJ,LJnodes,LJmode); % factor of two so that rc is node radius
        end
    end
end
forceArray = forceArray + Fc;