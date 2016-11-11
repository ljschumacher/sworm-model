function [ reflected ] = reflectDirection2D(angles,dimension)
% reflects movement direction upon incidence onto noflux boundary
% angles have to be given in [theta, phi], where theta and phi are column
% vectors
x =     1;
y =     2;

switch dimension
    case x % phi goes to pi - phi (or -pi + phi)
        angles(:,1) = wrapToPi(pi - angles(:,1));
    case y % phi goes to -phi
        angles(:,1) = - angles(:,1);
end
reflected = angles;
end

