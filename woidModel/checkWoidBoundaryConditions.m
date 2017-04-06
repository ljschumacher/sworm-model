function [ xyiarray ] = checkWoidBoundaryConditions( xyiarray, bc, L)
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
% phi =   3;

N = size(xyiarray,1); % number of objects
M = size(xyiarray,2); % number of nodes
ndim = size(xyiarray,3); % number of dims (x,y)

if iscell(bc)&&numel(bc)==ndim
    for dimCtr = [x y]
        switch bc{dimCtr}
            case 'periodic'
                nodeIndsUnder0 = find(xyiarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    xyiarray(nodeIndsUnder0)  = xyiarray(nodeIndsUnder0) + L(dimCtr);
                    nodeIndsOverL = find(xyiarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    xyiarray(nodeIndsOverL)  = xyiarray(nodeIndsOverL) - L(dimCtr);
                else % scalar domain size
                    xyiarray(nodeIndsUnder0)  = xyiarray(nodeIndsUnder0) + L;
                    nodeIndsOverL = find(xyiarray(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    xyiarray(nodeIndsOverL)  = xyiarray(nodeIndsOverL) - L;
                end
            case 'noflux'
                nodeIndsUnder0 = find(xyiarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                xyiarray(nodeIndsUnder0)  = - xyiarray(nodeIndsUnder0);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    nodeIndsOverL = find(xyiarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    xyiarray(nodeIndsOverL)  = 2*L(dimCtr) - xyiarray(nodeIndsOverL);
                else % scalar domain size
                    nodeIndsOverL = find(xyiarray(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    xyiarray(nodeIndsOverL)  = 2*L - xyiarray(nodeIndsOverL);
                end
%                 % change direction of movement upon reflection
%                 xyiarray(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr)) = ... % ugly use of indexing
%                     reflectDirection2D(...
%                     xyiarray(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr))...
%                     ,dimCtr);
        end
    end
else
    switch bc
        case 'periodic'
            if numel(L)==ndim % vector domain size [L_x L_y]
                for dimCtr = [x y]
                    nodeIndsUnder0 = find(xyiarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                    xyiarray(nodeIndsUnder0)  = xyiarray(nodeIndsUnder0) + L(dimCtr);
                    nodeIndsOverL = find(xyiarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    xyiarray(nodeIndsOverL)  = xyiarray(nodeIndsOverL) - L(dimCtr);
                end
            else % scalar domain size WARNING: periodic boundaries do not yet work with circular boundary
                nodeLogIndUnder0 = xyiarray(:,:,[x y])<0;
                xyiarray(nodeLogIndUnder0)  = xyiarray(nodeLogIndUnder0) + L;
                nodeLogIndOverL = xyiarray(:,:,[x y])>=L;
                xyiarray(nodeLogIndOverL)  = xyiarray(nodeLogIndOverL) - L;
            end
        case 'noflux'   
            if numel(L)==ndim % vector domain size [L_x L_y]
                for dimCtr = [x y]
                    nodeIndsUnder0 = find(xyiarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                    xyiarray(nodeIndsUnder0)  = - xyiarray(nodeIndsUnder0);
                    nodeIndsOverL = find(xyiarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    if any(nodeIndsOverL)
                        xyiarray(nodeIndsOverL)  = 2*L(dimCtr) - xyiarray(nodeIndsOverL);
                    end
%                     if any(nodeIndsUnder0)||any(nodeIndsOverL)
%                         % change direction of movement upon reflection
%                         xyiarray(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr)) = ... % ugly use of indexing
%                             reflectDirection2D(...
%                             xyiarray(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr))...
%                             ,dimCtr);
%                     end
                end
            else % scalar domain size --> circular domain boundary
                nodeIndsOverL = find(sqrt(sum(xyiarray(:,:,[x y]).^2,3))>=L);
                if any(nodeIndsOverL)
                    angles = atan2(xyiarray(nodeIndsOverL + N*M), ...%y coords
                    xyiarray(nodeIndsOverL)); %x coords
                    radii = sqrt(xyiarray(nodeIndsOverL + N*M).^2 + ...
                        xyiarray(nodeIndsOverL).^2);
                    xyiarray(nodeIndsOverL)  = (2*L - radii).*cos(angles);
                    xyiarray(nodeIndsOverL + N*M)  = (2*L - radii).*sin(angles);
%                     % change direction of movement upon reflection
%                     xyiarray(nodeIndsOverL + N*M*(phi - 1)) = ...
%                         alignWithBoundaryCircular(xyiarray(nodeIndsOverL + N*M*(phi - 1)),...
%                     angles);
                end
            end    
    end
end

end

