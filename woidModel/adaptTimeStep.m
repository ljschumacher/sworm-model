function [ dT ] = adaptTimeStep( dT0,v0,forceArray )
% adapt time-step such that it scales inversily with the max force
%   dT0 is a sensibly chosen starting time-step, e.g. rc/v0/8
%   v0 is target speed, and magnitude of motile forces, so that dT = dT0 if
%       motile forces are acting with intended magnitude

dT = dT0*v0/max(max(sqrt(sum(forceArray.^2,3))));

end

