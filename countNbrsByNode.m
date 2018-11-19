function [ numNbrsAndNodes ] = countNbrsByNode(distanceMatrix,r,nodeIndcs)
% counts how many nodes of each object have neighbours at r_ij<r
% distance matrix should have the form N by M by N by M with scalar
% distance values
if nargin<3
    nodeIndcs = 1:size(distanceMatrix,2);
end
N = size(distanceMatrix,1);
% check distance from each specified node (of each object) to any other nodes of all other objects
nbrNMNMLogInd = distanceMatrix(:,nodeIndcs,:,:)<=r;
for n = 1:N
    nbrNMNMLogInd(n,:,n,:) = false; % exclude self from neighbors
end
% nbrNMLogInd = any(nbrNMNLogInd,3); % don't care how many other objects in contact with
% numNodeswNbrs = sum(nbrNMLogInd,2); % count how many nodes of each object have neighbours
numNbrsPerNode = sum(sum(nbrNMNMLogInd,3),4); % count how many nodes of other objects each node is in contact with
numNbrsAndNodes = sum(numNbrsPerNode,2); % sum num nbrs over nodes (of this object)
end

