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
movState = 1 - 2*reversals(:,end); % =-1 if worm is reversing, 1 if not
Omega = (omega.*movState)*ones(1,M); % signed angular velocities for each worm and its nodes
reversalChanges = diff(reversals,1,2);
% if a reversal starts or ends, reset the phase based on current slope of
% shape at head
if any(reversalChanges ~=0)
    thetaNormalised = unwrap(thetaPrev,[],2) - mean(unwrap(thetaPrev,[],2),2); % subtract overall orientation
    dThetads = gradient(thetaNormalised,-1); % implicitly using |ds| = 1
    dThetadsNormalised = dThetads - mean(dThetads,2); % center on 0
    dThetadsNormalised = dThetadsNormalised./max(abs(dThetadsNormalised),[],2); % normalise range
    % use atan2 to get 4-quadrant angle back
    phaseFromShape = wrapTo2Pi(atan2(thetaNormalised,dThetadsNormalised)); % phase is given by atan(theta,dTheta/ds*l/deltaPhase)
    headIndcs = ~reversals(reversalChanges ~=0,end) + M*reversals(reversalChanges ~=0,end);
    allIndcsOrdered = ~reversals(reversalChanges ~=0,end)*(1:M) + reversals(reversalChanges ~=0,end)*(M:-1:1);
    try
    phaseOffset(reversalChanges ~=0,:) = wrapTo2Pi(phaseFromShape(reversalChanges ~=0,headIndcs)...
       - deltaPhase*(allIndcsOrdered - 1));
    catch
        error('phase reset went wrong')
    end
end
% 2nd order Runge-Kutta method with step-size h=1 (using gradient at
% midpoint to update)
thetaNow = wrapToPi(thetaPrev + theta_0*Omega.*cos(Omega*1/2 + phaseOffset)...
            + pi*reversalChanges); % 180 degree turn when reversal starts or ends
phaseOffset = wrapTo2Pi(phaseOffset + Omega*1); % update internal oscillator time by one time-step

end

