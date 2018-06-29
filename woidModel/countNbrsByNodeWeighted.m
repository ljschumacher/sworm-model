function [ numNbrsAndNodes ] = countNbrsByNodeWeighted(distanceMatrix,ri,rcontact,nodeIndcs)
% counts how many nodes of each object have neighbours at r_ij<r
% distance matrix should have the form N by M by N by M with scalar
% distance values
if nargin<3
    nodeIndcs = 1:size(distanceMatrix,2);
end
N = size(distanceMatrix,1);
distances = distanceMatrix(:,nodeIndcs,:,:);
% check distance from each specified node (of each object) to any other nodes of all other objects
nbrNMNMLogInd = distances<=ri;
for n = 1:N
    nbrNMNMLogInd(n,:,n,:) = false; % exclude self from neighbors
end
% weight neighbours by distance, such that any neighbours further than the
% contact distance get a 1/r weighting
nbrNMNMweighted = nbrNMNMLogInd.*rcontact./(max(distances,3*rcontact) - 2*rcontact);
numNbrsPerNode = sum(sum(nbrNMNMweighted,3),4); % count how many nodes of other objects each node is in contact with
numNbrsAndNodes = sum(numNbrsPerNode,2); % sum num nbrs over nodes (of this object)
end

