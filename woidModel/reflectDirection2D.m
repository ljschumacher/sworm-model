function [ reflected ] = reflectDirection2D(angles,dimension)
% reflects movement direction upon incidence onto noflux boundary

x =     1;
y =     2;

switch dimension
    case x % phi goes to pi - phi (or -pi + phi)
        reflected = wrapToPi(pi - angles);
    case y % phi goes to -phi
        reflected = - angles;
end
end

