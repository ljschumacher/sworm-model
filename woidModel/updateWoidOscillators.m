function thetaNow = updateWoidOscillators(thetaPrev, theta_0, omega_m, t, phaseOffset, reversals)
% updates the internal osciallators of woids
% INPUTS:
% thetaPrev - N by M matrix of internal oscillators at previous time step
% theta_0 - amplitude of oscillations
% omega_m - scalar or N by 1 vector of angular velocities
% t - scalar time point
% phaseOffset - scalar phase difference between adjacent nodes
% reverals - N by 1 logical index of which worms are reversing

M = size(thetaPrev,2);
movState = 1 - 2*reversals; % =-1 if worm is reversing, 1 if not
Omega = (omega_m.*movState)*ones(1,M); % signed angular velocities for each worm and its nodes
% 2nd order Runge-Kutta method with step-size h=1
thetaNow = thetaPrev + theta_0*Omega.*cos(Omega*(t-1/2) + phaseOffset);

end

