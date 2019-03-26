%Function for computing the euclidean distances between points, accounting
%for periodic boundary conditions

function D2 = periodiceucdist(XI,XJ)  

% Set the length of the period (applies to x and y directions)
period = 7.5;

dx = abs(XI-XJ);
per_dx = period-dx;
shortest = min(dx, per_dx);

sqdx = shortest.^2;
D2 = sqrt(sum(sqdx,2)); 
end

