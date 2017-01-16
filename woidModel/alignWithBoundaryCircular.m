function [ aligned ] = alignWithBoundaryCircular(angles,offSet)
% reflects movement direction upon incidence onto noflux boundary

aligned = wrapToPi(sign(angles - offSet)*pi/2 + offSet);

end

