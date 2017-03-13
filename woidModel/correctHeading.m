function thetaCorrected = correctHeading(xyarray,theta,v)

% corrects the heading of the worm based on how the worm actually moved
% if the actual displacement was in the same direction as the previous
% heading, the heading is not corrected
% if the displacement was different (for example because of collisions or
% other forces), the heading is retrospectively corrected
actualDisplacement = diff(xyarray,1,4);
targetDisplacement = v.*cat(3,cos(theta),sin(theta));
combined = actualDisplacement + targetDisplacement;
thetaCorrected = atan2(combined(:,:,2),combined(:,:,1));

end

