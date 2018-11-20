function headingsCorrected = correctHeading(forcearray)
% corrects the heading of the worm based on how the worm actually moved
% if the actual displacement was in the same direction as the previous
% heading, the heading is not corrected
% if the displacement was different (for example because of collisions or
% other forces), the heading is retrospectively corrected

% issues/to-do:
% - could be implemented in calculateForces instead, as effect of forces on
% orientation?

headingsCorrected = atan2(forcearray(:,:,2),forcearray(:,:,1));
end

