function thetaCorrected = correctHeading(xyarray,theta,v,bc,L)
% corrects the heading of the worm based on how the worm actually moved
% if the actual displacement was in the same direction as the previous
% heading, the heading is not corrected
% if the displacement was different (for example because of collisions or
% other forces), the heading is retrospectively corrected

% issues/to-do:
% - creates noisy headings at very slow speeds
% - does not work with mixed boundary conditions
% - could be implemented in calculateForces instead, as effect of forces on
% orientation

actualDisplacement = diff(xyarray,1,4);
if strcmp(bc,'periodic')
    for dimCtr = 1:length(L)
        actualDisplacement(:,:,dimCtr) = mod(actualDisplacement(:,:,dimCtr),L(dimCtr));
    end
end
targetDisplacement = v.*cat(3,cos(theta),sin(theta));
combined = actualDisplacement + targetDisplacement;
thetaCorrected = atan2(combined(:,:,2),combined(:,:,1));

end

