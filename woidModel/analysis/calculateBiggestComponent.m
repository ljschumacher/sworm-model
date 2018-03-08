function [biggestComponentSize] = calculateBiggestComponent(positions,interactionRadius)
% calculates the biggest component in simulations, based on objects in
% contact
% INPUTS
% positions is of size N,M,2

distanceMatrix = sqrt(sum(computeWoidDistancesWithBCs(positions,[],'free').^2,5));
adjacencyMatrix = squeeze(any(any(distanceMatrix<=interactionRadius,4),2)); % worms 'touching' at any of their nodes
connectedComponents = conncomp(graph(adjacencyMatrix));
biggestComponentSize = nnz(connectedComponents==mode(connectedComponents));
                
end

