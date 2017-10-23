function [ v, omega ] = slowWorms(distanceMatrix,ri,slowingNodes,slowingMode,...
    vs,v0,omega_m,num_nbr_max_per_node,roamingLogInd)
% slow down worms based on contact with other worms
% slowLogInd = findWoidNeighbors(distanceMatrix,ri,slowingNodes);
% v = vs*slowLogInd + v0*(~slowLogInd); % adjust speed for slowed worms
% omega = omega_m*(~slowLogInd + vs/v0*slowLogInd); % adjust internal oscillator freq for slowed worms

if ~isempty(slowingNodes)&&any(~roamingLogInd)
    if strcmp(slowingMode,'density')
        num_nbr_max = numel(slowingNodes)*num_nbr_max_per_node;
        num_nbr = countWoidNeighbors(distanceMatrix,ri,slowingNodes); % number of nodes in contact with neighbours
        num_nbr = min(num_nbr,num_nbr_max); % depending on how we count neighbours, we may have more than the max number
        % slow worms proportional to how many nbrs it has cumulatively over all nodes
        v = v0*(1 - num_nbr/num_nbr_max) + vs*(num_nbr/num_nbr_max);
    else
        m_nbr_max = numel(slowingNodes);
        m_nbr = countNodesContacted(distanceMatrix,ri,slowingNodes); % number of nodes in contact with neighbours
        if strcmp(slowingMode,'gradual')
            % slow worms proportional to how many nodes are in contact with other nodes
            v = v0*(1 - m_nbr/m_nbr_max) + vs*(m_nbr/m_nbr_max);
        else
            slowWorms = m_nbr>0;
            v = v0*~slowWorms + vs*slowWorms;
        end
    end
    v(roamingLogInd) = v0; % roaming worms don't slow down
    assert(~any(v>v0+eps(v0)))
    omega = omega_m*v/v0; % adjust internal oscillator freq for slowed worms
else
    N = size(distanceMatrix,1);
    v = v0*ones(N,1);
    omega = omega_m;
end
end

