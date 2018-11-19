function [ aligned ] = alignWithBoundary2D(heading,dimension)
% reflects movement direction upon incidence onto noflux boundary

x =     1;
y =     2;

switch dimension
    case x % heading goes to +-pi/2
        aligned = sign(wrapToPi(2*heading)/2)*pi/2;
    case y % heading goes to +-pi
        aligned = sign(abs(heading) - pi/2)*pi;
end
end

