function [ aligned ] = alignWithBoundary2D(angles,dimension)
% reflects movement direction upon incidence onto noflux boundary

x =     1;
y =     2;

switch dimension
    case x % phi goes to +-pi/2
        aligned = sign(wrapToPi(2*angles)/2)*pi/2;
    case y % phi goes to +-pi
        aligned = sign(abs(angles) - pi/2)*pi;
end
end

