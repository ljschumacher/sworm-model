function [ xyarray, theta ] = checkWoidBoundaryConditions(xyarray, theta, bc, L)
% check boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can
%   be single number or 2 element array {'bcx','bcy'} for different
%   bcs along different dimensions

% issues/to-do's:
%   - could try to optimise the code for the vector domain size/mixed
%   boundary condition cases to get rid of loops...
%   - mixed periodic boundary conditions can be quite slow?

% short-hand for indexing coordinates
x =     1;
y =     2;

N = size(xyarray,1); % number of objects
M = size(xyarray,2); % number of nodes
ndim = size(xyarray,3); % number of dims (x,y)

if iscell(bc)&&numel(bc)==ndim
    for dimCtr = [x y]
        switch bc{dimCtr}
            case 'periodic'
                nodeIndsUnder0 = find(xyarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    xyarray(nodeIndsUnder0)  = mod(xyarray(nodeIndsUnder0),L(dimCtr));
                    nodeIndsOverL = find(xyarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    xyarray(nodeIndsOverL)  = mod(xyarray(nodeIndsOverL),L(dimCtr));
                else % scalar domain size
                    xyarray(nodeIndsUnder0)  = mod(xyarray(nodeIndsUnder0),L);
                    nodeIndsOverL = find(xyarray(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    xyarray(nodeIndsOverL)  = mod(xyarray(nodeIndsOverL),L);
                end
            case 'noflux'
                nodeIndsUnder0 = find(xyarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                xyarray(nodeIndsUnder0)  = - xyarray(nodeIndsUnder0);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    nodeIndsOverL = find(xyarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    xyarray(nodeIndsOverL)  = 2*L(dimCtr) - xyarray(nodeIndsOverL);
                else % scalar domain size
                    nodeIndsOverL = find(xyarray(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    xyarray(nodeIndsOverL)  = 2*L - xyarray(nodeIndsOverL);
                end
                % change direction of movement upon reflection
                theta(union(nodeIndsUnder0,nodeIndsOverL) - N*M*(dimCtr - 1)) = ... % ugly use of indexing
                    reflectDirection2D(theta(union(nodeIndsUnder0,nodeIndsOverL) - N*M*(dimCtr - 1)),dimCtr);
        end
    end
else
    switch bc
        case 'periodic'
            if numel(L)==ndim % vector domain size [L_x L_y]
                for dimCtr = [x y]
                    nodeIndsUnder0 = find(xyarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1); % don't use logical indexing as we want to allow for non-square domains
                    xyarray(nodeIndsUnder0) = mod(xyarray(nodeIndsUnder0),L(dimCtr));
                    nodeIndsOverL = find(xyarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    xyarray(nodeIndsOverL) = mod(xyarray(nodeIndsOverL),L(dimCtr));
                end
            else % scalar domain size WARNING: periodic boundaries do not yet work with circular boundary
                nodeLogIndUnder0 = xyarray(:,:,[x y])<0;
                nodeLogIndOverL = xyarray(:,:,[x y])>=L;
                xyarray(nodeLogIndUnder0|nodeLogIndOverL)  = mod(xyarray(nodeLogIndUnder0),L);
            end
        case 'noflux'
            if numel(L)==ndim % vector domain size [L_x L_y]
                for dimCtr = [x y]
                    nodeIndsUnder0 = find(xyarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                    xyarray(nodeIndsUnder0)  = - xyarray(nodeIndsUnder0);
                    nodeIndsOverL = find(xyarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    if any(nodeIndsOverL)
                        xyarray(nodeIndsOverL)  = 2*L(dimCtr) - xyarray(nodeIndsOverL);
                    end
                    if any(nodeIndsUnder0)||any(nodeIndsOverL)
                        % change direction of movement upon reflection
                        theta(union(nodeIndsUnder0,nodeIndsOverL) - N*M*(dimCtr - 1)) = ... % ugly use of indexing
                            reflectDirection2D(theta(union(nodeIndsUnder0,nodeIndsOverL) - N*M*(dimCtr - 1)),dimCtr);
                    end
                end
            else % scalar domain size --> circular domain boundary
                nodeIndsOverL = find(sqrt(sum(xyarray(:,:,[x y]).^2,3))>=L);
                if any(nodeIndsOverL)
                    angles = atan2(xyarray(nodeIndsOverL + N*M), ...%y coords
                        xyarray(nodeIndsOverL)); %x coords
                    radii = sqrt(xyarray(nodeIndsOverL + N*M).^2 + ...
                        xyarray(nodeIndsOverL).^2);
                    xyarray(nodeIndsOverL)  = (2*L - radii).*cos(angles);
                    xyarray(nodeIndsOverL + N*M)  = (2*L - radii).*sin(angles);
%                     % change direction of movement upon reflection
%                     theta(nodeIndsOverL) = ...
%                         alignWithBoundaryCircular(theta(nodeIndsOverL),...
%                         angles);
                end
            end
    end
end

end

