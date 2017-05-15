function [ reflected ] = reflectDirection2D(heading,dimension)
% reflects movement direction upon incidence onto noflux boundary

x =     1;
y =     2;

switch dimension
    case x % heading goes to pi - heading (or -pi + heading)
        reflected = wrapToPi(pi - heading);
    case y % heading goes to -heading
        reflected = - heading;
end
end

