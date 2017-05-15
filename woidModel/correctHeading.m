function headingsCorrected = correctHeading(forcearray,headings,bc,L)
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
% N = size(xyarray,1);
% M = size(xyarray,2);
% 
% actualDisplacement = diff(xyarray,1,4);
% if strcmp(bc,'periodic')
%     for dimCtr = 1:length(L)
%         overIndcs = find(actualDisplacement(:,:,dimCtr)>=L(dimCtr)/2) + N*M*(dimCtr - 1);
%         actualDisplacement(overIndcs) = actualDisplacement(overIndcs) - L(dimCtr);
%         underIndcs = find(actualDisplacement(:,:,dimCtr)<=-L(dimCtr)/2) + N*M*(dimCtr - 1);
%         actualDisplacement(underIndcs) = actualDisplacement(underIndcs) + L(dimCtr);
%     end
% end
% targetDisplacement = v.*cat(3,cos(headings),sin(headings));
% combined = actualDisplacement;% + targetDisplacement;
% headingsCorrected = atan2(combined(:,:,2),combined(:,:,1));
headingsCorrected = atan2(forcearray(:,:,2),forcearray(:,:,1));
end

