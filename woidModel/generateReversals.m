function reversalLogInd = generateReversals(reversalLogInd,timeCtr,distanceMatrix,...
    interactionRadius,contactRadius,headInd,tailInd,dT,reversalMode,revRate,revTime,revRateCluster,...
    revRateClusterEdge,roamingLogInd,drdN_rev)
if strcmp(reversalMode,'density')
    % increase reversal rate based on local density
    num_nbr_head = countNbrsByNode(distanceMatrix,interactionRadius,headInd)./numel(headInd); % average number of neighbouring nodes in contact
    num_nbr_tail = countNbrsByNode(distanceMatrix,interactionRadius,tailInd)./numel(tailInd); % average number of neighbouring nodes in contact
    revRateClusterEdgeTailOut = revRateClusterEdge + drdN_rev.*num_nbr_head;
    revRateClusterEdgeHeadOut = revRateClusterEdge + drdN_rev.*num_nbr_tail;
elseif strcmp(reversalMode,'density_weighted')
    % increase reversal rate based on local density
    num_nbr_head = countNbrsByNodeWeighted(distanceMatrix,interactionRadius,contactRadius,headInd)./numel(headInd); % average number of neighbouring nodes in contact
    num_nbr_tail = countNbrsByNodeWeighted(distanceMatrix,interactionRadius,contactRadius,tailInd)./numel(tailInd); % average number of neighbouring nodes in contact
    revRateClusterEdgeTailOut = revRateClusterEdge + drdN_rev.*num_nbr_head;
    revRateClusterEdgeHeadOut = revRateClusterEdge + drdN_rev.*num_nbr_tail;
end

% find which worms are poking their head or tail out
tailContacts = findWoidNeighbors(distanceMatrix,3/2*contactRadius,tailInd);
headContacts = findWoidNeighbors(distanceMatrix,3/2*contactRadius,headInd);

% generates reversal states
% find worms currently in reversal state
currentReversalsLogInd = reversalLogInd(:,timeCtr);

% worms outside clusters reverse with revRate for revTime, unless already reversing
freeFwdWormsLogInd = ~tailContacts&~headContacts&~currentReversalsLogInd;
reversalLogInd(freeFwdWormsLogInd,timeCtr:(timeCtr+revTime)) ... % set reversal state for duration of reversal
    = repmat(logical(poissrnd(revRate*dT,nnz(freeFwdWormsLogInd),1)),1,revTime+1);

% worms inside clusters reverse with revRateCluster for revTime, unless already reversing
clustFwdWormsLogInd = tailContacts&headContacts&~currentReversalsLogInd;
reversalLogInd(clustFwdWormsLogInd,timeCtr:(timeCtr+revTime)) ... % set reversal state for duration of reversal
    = repmat(logical(poissrnd(revRateCluster*dT,nnz(clustFwdWormsLogInd),1)),1,revTime+1);

% stop reversal with increased rate if tail is sticking out of cluster (unless roaming)
freeBwdTailsInd = find(~tailContacts&headContacts&currentReversalsLogInd&~roamingLogInd); % Ntrue by 1
if strcmp(reversalMode,'contact')
    stoppedReversalsLogInd = logical(poissrnd(revRateClusterEdge*dT,numel(freeBwdTailsInd),1)); % Ntrue by 1
elseif strcmp(reversalMode,'density')||strcmp(reversalMode,'density_weighted')
    stoppedReversalsLogInd = logical(poissrnd(revRateClusterEdgeTailOut(freeBwdTailsInd)*dT,numel(freeBwdTailsInd),1)); % Ntrue by 1
end
reversalLogInd(freeBwdTailsInd(stoppedReversalsLogInd),timeCtr:end) = false;

% reverse with increased rate if head is sticking out of cluster (unless
% roaming)
freeFwdHeadsLogInd = tailContacts&~headContacts&~currentReversalsLogInd&~roamingLogInd; % N by 1
if strcmp(reversalMode,'contact')
    reversalLogInd(freeFwdHeadsLogInd,timeCtr:(timeCtr+revTime)) ...
        = repmat(logical(poissrnd(revRateClusterEdge*dT,nnz(freeFwdHeadsLogInd),1)),1,revTime+1);
elseif strcmp(reversalMode,'density')||strcmp(reversalMode,'density_weighted')
    reversalLogInd(freeFwdHeadsLogInd,timeCtr:(timeCtr+revTime)) ...
        = repmat(logical(poissrnd(revRateClusterEdgeHeadOut(freeFwdHeadsLogInd)*dT,nnz(freeFwdHeadsLogInd),1)),1,revTime+1);
end
end

