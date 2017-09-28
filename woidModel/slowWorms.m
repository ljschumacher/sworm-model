function [ v, omega ] = slowWorms(distanceMatrix,ri,slowingNodes,slowingMode,vs,v0,omega_m)
% slow down worms based on contact with other worms
% slowLogInd = findWoidNeighbors(distanceMatrix,ri,slowingNodes);
% v = vs*slowLogInd + v0*(~slowLogInd); % adjust speed for slowed worms
% omega = omega_m*(~slowLogInd + vs/v0*slowLogInd); % adjust internal oscillator freq for slowed worms

m_nbr_max = numel(slowingNodes);
if m_nbr_max > 0
    m_nbr = countWoidNeighbors(distanceMatrix,ri,slowingNodes); % number of nodes in contact with neighbours
    if strcmp(slowingMode,'gradual')
        % slow worms proportional to how many nodes are in contact with other nodes
    v = v0*(1 - m_nbr/m_nbr_max) + vs*(m_nbr/m_nbr_max);
    else
        slowWorms = m_nbr>0;
        v = v0*~slowWorms + vs*slowWorms;
    end
    assert(~any(v>v0+eps(v0)))
    omega = omega_m*v/v0; % adjust internal oscillator freq for slowed worms
else
    N = size(distanceMatrix,1);
    v = v0*ones(N,1);
    omega = omega_m;
end
end

