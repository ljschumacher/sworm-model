function [ v, omega ] = slowWorms(distanceMatrix,ri,slowingNodes,vs,v0,omega_m)
% slow down worms based on contact with other worms
slowLogInd = findWoidNeighbors(distanceMatrix,2*ri,slowingNodes);
v = vs*slowLogInd + v0*(~slowLogInd); % adjust speed for slowed worms
omega = omega_m*(~slowLogInd + vs/v0*slowLogInd); % adjust internal oscillator freq for slowed worms
end

