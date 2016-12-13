function reversalLogInd = generateReversals(reversalLogInd,t,distanceMatrix,...
    contactRadius,headInd,tailInd,revRate,revTime,revRateReduced,revTimeReduced)
% find which worms are in or out of clusters    
tailContacts = findWoidNeighbors(distanceMatrix,contactRadius,tailInd);
headContacts = findWoidNeighbors(distanceMatrix,contactRadius,headInd);
% generates reversal states
% find worms currently in reversal state
currentReversals = reversalLogInd(:,t);
% worms outside clusters reverse with revRate for revTime, unless already reversing
freeFwdWorms = ~tailContacts&~headContacts&~currentReversals;
reversalLogInd(freeFwdWorms,t:(t+revTime)) ... % set reversal state for duration of reversal
    = repmat(logical(poissrnd(revRate,nnz(freeFwdWorms),1)),1,revTime+1);
% worms inside clusters reverse with revRateReduced for revTimeReduced, unless already reversing
clustFwdWorms = tailContacts&headContacts&~currentReversals;
reversalLogInd(clustFwdWorms,t:(t+revTimeReduced)) ... % set reversal state for duration of reversal
    = repmat(logical(poissrnd(revRateReduced,nnz(clustFwdWorms),1)),1,revTimeReduced+1);

% stop reversal if tail is sticking out of cluster
freeBwdTails = ~tailContacts&headContacts&currentReversals;
reversalLogInd(freeBwdTails,t:end) = false;
% reverse with normal rate if head is sticking out of cluster
freeFwdHeads = tailContacts&~headContacts&~currentReversals;
reversalLogInd(freeFwdHeads,t:(t+revTime)) ...
    = repmat(logical(poissrnd(revRate,nnz(freeFwdHeads),1)),1,revTime+1);

end

