function [thetaNow, phaseOffset] = updateWoidOscillators(thetaPrev, theta_0, ...
    omega, phaseOffset, deltaPhase, reversals)
% updates the internal osciallators of woids
% INPUTS:
% thetaPrev - N by M matrix of internal oscillators at previous time step
% theta_0 - amplitude of oscillations
% omega_m - scalar or N by 1 vector of angular velocities
% t - scalar time point
% phaseOffset - scalar phase difference between adjacent nodes
% reverals - N by 1 logical index of which worms are reversing

M = size(thetaPrev,2);
movState = 1 - 2*reversals(:,2); % =-1 if worm is reversing, 1 if not
Omega = (omega.*movState)*ones(1,M); % signed angular velocities for each worm and its nodes
reversalChanges = diff(reversals,1,2);
% 2nd order Runge-Kutta method with step-size h=1 (using gradient at
% midpoint to update)
thetaNow = wrapToPi(thetaPrev + theta_0*Omega.*cos(Omega*1/2 + phaseOffset)...
            + pi*reversalChanges); % 180 degree turn when reversal starts or ends
phaseOffset = wrapTo2Pi(phaseOffset + Omega*1); % update internal oscillator time by one time-step
% if a reversal starts or ends, reset the phase based on current slope of
% shape at head
% if any(reversalChanges ~=0)
%     dThetads = gradient(unwrap(thetaNow,[],2)); % implicitly using |ds| = 1
%     dThetadsNormalised = dThetads - mean(dThetads,2);
%     phaseFromShape = wrapTo2Pi(acos(dThetadsNormalised./max(abs(dThetadsNormalised),[],2)));
%     headIndcs = ~reversals(reversalChanges ~=0) + M*reversals(reversalChanges ~=0);
%     allIndcsOrdered = ~reversals(reversalChanges ~=0)*(1:M) + reversals(reversalChanges ~=0)*(M:-1:1);
%     phaseOffset(reversalChanges ~=0,:) = wrapTo2Pi(phaseFromShape(reversalChanges ~=0,headIndcs)...
%         + deltaPhase*allIndcsOrdered);
% end
% % evaluate gradient along arc length for t+1/2
% dThetadsPrev = gradient(unwrap(thetaPrev));
% thetaMidpoint = wrapToPi(thetaPrev +  Omega./deltaPhase.*dThetadsPrev);
% dThetadsMidpoint = gradient(unwrap(thetaMidpoint));
% thetaNow = wrapToPi(thetaPrev +  Omega./deltaPhase.*dThetadsMidpoint);
end

