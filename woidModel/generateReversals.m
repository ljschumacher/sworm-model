function reversalLogInd = generateReversals(reversalLogInd,timeCtr,distanceMatrix,...
    interactionRadius,headInd,tailInd,dT,revRate,revTime,revRateCluster,revRateClusterEdge,roamingLogInd)
% find which worms are in or out of clusters    
tailContacts = findWoidNeighbors(distanceMatrix,interactionRadius,tailInd);
headContacts = findWoidNeighbors(distanceMatrix,interactionRadius,headInd);
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

% stop reversal with increased rate if tail is sticking out of cluster
% (unless roaming)
freeBwdTailsInd = find(~tailContacts&headContacts&currentReversalsLogInd&~roamingLogInd); % Ntrue by 1
stoppedReversalsLogInd = logical(poissrnd(revRateClusterEdge*dT,numel(freeBwdTailsInd),1)); % Ntrue by 1
reversalLogInd(freeBwdTailsInd(stoppedReversalsLogInd),timeCtr:end) = false;
% reverse with increased rate if head is sticking out of cluster (unless
% roaming)
freeFwdHeadsLogInd = tailContacts&~headContacts&~currentReversalsLogInd&~roamingLogInd; % N by 1
reversalLogInd(freeFwdHeadsLogInd,timeCtr:(timeCtr+revTime)) ...
    = repmat(logical(poissrnd(revRateClusterEdge*dT,nnz(freeFwdHeadsLogInd),1)),1,revTime+1);

end

