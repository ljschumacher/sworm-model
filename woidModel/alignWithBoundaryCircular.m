function [ aligned ] = alignWithBoundaryCircular(heading,offSet)
% reflects movement direction upon incidence onto noflux boundary

aligned = wrapToPi(sign(heading - offSet)*pi/2 + offSet);

end

