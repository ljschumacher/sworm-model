function [positions, anglesInitial] = initialiseWoids(N,M,numTimepoints,L,segmentLength,phaseOffset,theta_0,exclusionRadius,bc)
% initialises object positions and directions
% uniformly randomly distributed

% issues/to-do's:
% - needlessy re-computing distances between all successfully initialised
% agents when checking for overlaps
% - not tested for periodic boundaries

positions = NaN(N,M,2,numTimepoints);
anglesInitial = NaN(N,M);
% short-hand for indexing coordinates
x =     1;
y =     2;

if size(L,1)>size(L,2)
    L = L';
end

L0 = (M-1)*segmentLength;
disp(['Initialising woid positions...' datestr(now)])
objCtr = 1;
tryCtr = 1;
while objCtr <= N
    % initialise head node
    if numel(L)==1
        % Position, within circle of radius L-L0
        rndangle = pi*(2*rand(1) - 1);
        rndradius = (L - L0).*rand(1);
        positions(objCtr,1,[x y],1) = rndradius.*[cos(rndangle) sin(rndangle)]; % initialise positions at least one woid length away from edge
    else
        % Position within rectangle of [Lx, Ly]
        positions(objCtr,1,[x y],1) = L0 + (L - 2*L0).*rand(size(L));
    end
    % Direction
    anglesInitial(objCtr,:) = wrapToPi(pi*(2*rand - 1) + ...
        theta_0*sin(phaseOffset(objCtr,:)));   % random orientation between -pi and pi for each object plus undulations
    for nodeCtr = 2:M % initialise woid positions, node by node
        % initialise the next node in the right direction at segmentLength away
        positions(objCtr,nodeCtr,[x y],1) = squeeze(positions(objCtr,nodeCtr - 1,[x y],1))... % previous node's position
            - segmentLength*[cos(anglesInitial(objCtr,nodeCtr)); sin(anglesInitial(objCtr,nodeCtr))];
    end
    if exclusionRadius>0
        % check for any overlaps, and if so discard this woid's position
        distanceMatrix = sqrt(sum(computeWoidDistancesWithBCs(positions(1:objCtr,:,[x y],1),L,bc).^2,5));% scalar distance from current object to any other
        collisionNbrs = squeeze(any(any(2*exclusionRadius>=distanceMatrix(objCtr,:,:,:),4),2));
        collisionNbrs(objCtr) = false; % don't volume-exclude self
    else
        collisionNbrs=false;
    end
    if ~any(collisionNbrs)
        %         disp(['initialised ' num2str(objCtr) ' out of ' num2str(N) ' worms'...
        %             ' after ' num2str(tryCtr) ' tries'])
        objCtr = objCtr + 1;
        tryCtr = 1;
    else
        tryCtr = tryCtr + 1;
    end
end

end