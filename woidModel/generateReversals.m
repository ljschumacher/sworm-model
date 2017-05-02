function reversalLogInd = generateReversals(reversalLogInd,timeCtr,distanceMatrix,...
    interactionRadius,headInd,tailInd,dT,revRate,revTime,revRateReduced)
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
% worms inside clusters reverse with revRateReduced for revTime, unless already reversing
clustFwdWorms = tailContacts&headContacts&~currentReversals;
reversalLogInd(clustFwdWorms,timeCtr:(timeCtr+revTime)) ... % set reversal state for duration of reversal
    = repmat(logical(poissrnd(revRateReduced*dT,nnz(clustFwdWorms),1)),1,revTime+1);

% stop reversal if tail is sticking out of cluster
freeBwdTails = ~tailContacts&headContacts&currentReversals;
reversalLogInd(freeBwdTails,timeCtr:end) = false;
% reverse with normal rate if head is sticking out of cluster
freeFwdHeads = tailContacts&~headContacts&~currentReversals;
reversalLogInd(freeFwdHeads,timeCtr:(timeCtr+revTime)) ...
    = repmat(logical(poissrnd(revRate*dT,nnz(freeFwdHeads),1)),1,revTime+1);

end

