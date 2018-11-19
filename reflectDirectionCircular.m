function [ reflected ] = reflectDirectionCircular(angles,offSet)
% reflects movement direction upon incidence onto noflux boundary

reflected = wrapToPi(pi - angles + 2*offSet);

end

