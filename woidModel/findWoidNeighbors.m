function [ nbrLogInd ] = findWoidNeighbors(distanceMatrix,r)
% computes which object have neighbors at r_ij<r
% distance matrix should have the form N by M by N by M with scaler
% distance values

N = size(distanceMatrix,1);
% check distance from any nodes (of each object) to any other nodes of all other objects
nbrNNLogInd = squeeze(any(any(distanceMatrix<=r,4),2)); % reduce to N by N
nbrNNLogInd(logical(speye(N))) = false; % exclude self from neighbors
nbrLogInd = any(nbrNNLogInd,2); % find if each object has any neighbors
end

