function [ nbrNums ] = countWoidNeighbors(distanceMatrix,r,nodeIndcs)
% computes number of neighbours for each object at r_ij<r
% distance matrix should have the form N by M by N by M with scalar
% distance values
if nargin<3
    nodeIndcs = 1:size(distanceMatrix,2);
end
N = size(distanceMatrix,1);
% check distance from each node (of each object) to any other nodes of all other objects
nbrNMNLogInd = squeeze(any(distanceMatrix(:,nodeIndcs,:,:)<=r,4));
for n = 1:N
    nbrNMNLogInd(n,:,n) = false; % exclude self from neighbors
end
nbrNMLogInd = any(nbrNMNLogInd(n,:,n),3); 
nbrNums = sum(nbrNMLogInd,2); % count how many nodes of each object have neighbours
end

