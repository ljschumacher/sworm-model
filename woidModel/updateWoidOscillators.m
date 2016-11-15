function thetaNow = updateWoidOscillators(thetaPrev, theta_0, omega_m, t, phaseOffset, reversals)
% updates the internal osciallators of woids
M = size(thetaPrev,2);
movState = 1 - 2*reversals; % =-1 if worm is reversing, 1 if not
Omega = omega_m*movState*ones(1,M); % signed angular velocities for each worm and its nodes
thetaNow = thetaPrev + theta_0*Omega.*cos(Omega*t + phaseOffset);

end

