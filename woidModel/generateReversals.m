function reversalLogInd = generateReversals(reversalLogInd,timeCtr,distanceMatrix,...
    interactionRadius,headInd,tailInd,dT,revRate,revTime,revRateCluster,revRateClusterEdge,roamingLogInd)
% find which worms are in or out of clusters    
tailContacts = findWoidNeighbors(distanceMatrix,interactionRadius,tailInd);
headContacts = findWoidNeighbors(distanceMatrix,interactionRadius,headInd);
% generates reversal states
% find worms currently in reversal state
currentReversals = reversalLogInd(:,timeCtr);
% worms outside clusters reverse with revRate for revTime, unless already reversing
freeFwdWorms = ~tailContacts&~headContacts&~currentReversals;
reversalLogInd(freeFwdWorms,timeCtr:(timeCtr+revTime)) ... % set reversal state for duration of reversal
    = repmat(logical(poissrnd(revRate*dT,nnz(freeFwdWorms),1)),1,revTime+1);
% worms inside clusters reverse with revRateCluster for revTime, unless already reversing
clustFwdWorms = tailContacts&headContacts&~currentReversals;
reversalLogInd(clustFwdWorms,timeCtr:(timeCtr+revTime)) ... % set reversal state for duration of reversal
    = repmat(logical(poissrnd(revRateCluster*dT,nnz(clustFwdWorms),1)),1,revTime+1);

% stop reversal with increased rate if tail is sticking out of cluster
% (unless roaming)
freeBwdTails = ~tailContacts&headContacts&currentReversals&~roamingLogInd;
stoppedReversalsLogInd = logical(poissrnd(revRateClusterEdge*dT,nnz(freeBwdTails),1));
reversalLogInd(freeBwdTails(stoppedReversalsLogInd),timeCtr:end) = false;
% reverse with increased rate if head is sticking out of cluster (unless
% roaming)
freeFwdHeads = tailContacts&~headContacts&~currentReversals&~roamingLogInd;
reversalLogInd(freeFwdHeads,timeCtr:(timeCtr+revTime)) ...
    = repmat(logical(poissrnd(revRateClusterEdge*dT,nnz(freeFwdHeads),1)),1,revTime+1);

end

