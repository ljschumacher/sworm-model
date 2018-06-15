function [intersectCheck] = checkIntersects(positions,distanceMatrix,r_collision,bc,L)
% checks if any two worm chains cross, using the file exchange contribution
% interX

N = size(distanceMatrix,1);
intersectCheck = false;

% for efficieny we will only check chains that are in close proximity
collisionNbrsNMNM = distanceMatrix<r_collision; % check distance to all other nodes of all other objects
% collisionNbrsNMNM(objCtr,:) = false; % no contact force with self
if any(collisionNbrsNMNM(:)) % effectively only checks if rc>0, as it includes self-contacts within worm chain
    for ii = 1:N
        thisCollisionNbrs = any(collisionNbrsNMNM(:,:,ii,:),4); % take collisions from with any nodes of this worm
        thisCollisionNbrs = any(thisCollisionNbrs,2); % take collisions with any node of another worm
        thisCollisionNbrs(ii) = false; % do not check for intersects with self
        if any(thisCollisionNbrs(:))
            collNbrIdcs = find(thisCollisionNbrs)';
            collNbrIdcs2check = collNbrIdcs(collNbrIdcs>ii); % only check jj>ii, because of symmetry
            for jj=collNbrIdcs2check
                % need to undo periodic boundaries before checking for intersects
                if strcmp(bc,'periodic') % might fail for mixed boundary conditions
                    centeredPositions = centerPositions(positions([ii jj],:,:),L);
                    interxP = InterX(squeeze(centeredPositions(1,:,:))',squeeze(centeredPositions(2,:,:))'); % calculate intersection point(s) of two worms
                else
                    interxP = InterX(squeeze(positions(ii,:,:))',squeeze(positions(jj,:,:))'); % calculate intersection point(s) of two worms
                end
                if ~isempty(interxP)
                    intersectCheck = true;
                    break
                end
                
            end
        end
        if intersectCheck==true
            break
        end
    end
end

end

function positions = centerPositions(positions,L)
% this calculates the centre of mass for peridic boundaries - useful
% trick found on wikipedia
c_x = mean(mean(cos(positions(:,:,1)/L(1)*2*pi),2),1);
s_x = mean(mean(sin(positions(:,:,1)/L(1)*2*pi),2),1);
xoffset = L(2)/2/pi*(atan2(-s_x,-c_x) + pi) - L(1)/2;
c_y = mean(mean(cos(positions(:,:,2)/L(2)*2*pi),2),1);
s_y = mean(mean(sin(positions(:,:,2)/L(2)*2*pi),2),1);
yoffset = L(2)/2/pi*(atan2(-s_y,-c_y) + pi) - L(2)/2;
positions(:,:,1) = positions(:,:,1) - xoffset;
positions(:,:,2) = positions(:,:,2) - yoffset;
% re-enforce periodic boundaries
[ positions, ~ ] = checkWoidBoundaryConditions(positions, [], 'periodic', L);
end