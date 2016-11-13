function [ xyphiarray ] = checkWoidBoundaryConditions( xyphiarray, bc, L)
% check boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can
%   be single number or 2 element array {'bcx','bcy'} for different
%   bcs along different dimensions

% issues/to-do's:
%   - could try to optimise the code for the vector domain size/mixed
%   boundary condition cases to get rid of loops...
%   - mixed periodic boundary conditions can be quite slow?
%   - move check of boundary conditions to separate function?

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

N = size(xyphiarray,1); % number of objects
M = size(xyphiarray,2); % number of nodes
ndim = size(xyphiarray,3)-1; % number of dims (x,y)

if iscell(bc)&&numel(bc)==ndim
    for dimCtr = [x y]
        switch bc{dimCtr}
            case 'periodic'
                nodeIndsUnder0 = find(xyphiarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    xyphiarray(nodeIndsUnder0)  = xyphiarray(nodeIndsUnder0) + L(dimCtr);
                    nodeIndsOverL = find(xyphiarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    xyphiarray(nodeIndsOverL)  = xyphiarray(nodeIndsOverL) - L(dimCtr);
                else % scalar domain size
                    xyphiarray(nodeIndsUnder0)  = xyphiarray(nodeIndsUnder0) + L;
                    nodeIndsOverL = find(xyphiarray(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    xyphiarray(nodeIndsOverL)  = xyphiarray(nodeIndsOverL) - L;
                end
            case 'noflux'
                nodeIndsUnder0 = find(xyphiarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                xyphiarray(nodeIndsUnder0)  = - xyphiarray(nodeIndsUnder0);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    nodeIndsOverL = find(xyphiarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    xyphiarray(nodeIndsOverL)  = 2*L(dimCtr) - xyphiarray(nodeIndsOverL);
                else % scalar domain size
                    nodeIndsOverL = find(xyphiarray(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    xyphiarray(nodeIndsOverL)  = 2*L - xyphiarray(nodeIndsOverL);
                end
                % change direction of movement upon reflection
                xyphiarray(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr)) = ... % ugly use of indexing
                    reflectDirection2D(...
                    xyphiarray(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr))...
                    ,dimCtr);
        end
    end
else
    switch bc
        case 'periodic'
            if numel(L)==ndim % vector domain size [L_x L_y]
                for dimCtr = [x y]
                    nodeIndsUnder0 = find(xyphiarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                    xyphiarray(nodeIndsUnder0)  = xyphiarray(nodeIndsUnder0) + L(dimCtr);
                    nodeIndsOverL = find(xyphiarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    xyphiarray(nodeIndsOverL)  = xyphiarray(nodeIndsOverL) - L(dimCtr);
                end
            else % scalar domain size
                nodeLogiUnder0 = xyphiarray(:,:,[x y])<0;
                xyphiarray(nodeLogiUnder0)  = xyphiarray(nodeLogiUnder0) + L;
                nodeLogiOverL = xyphiarray(:,:,[x y])>=L;
                xyphiarray(nodeLogiOverL)  = xyphiarray(nodeLogiOverL) - L;
            end
        case 'noflux'
            for dimCtr = [x y]
                nodeIndsUnder0 = find(xyphiarray(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                xyphiarray(nodeIndsUnder0)  = - xyphiarray(nodeIndsUnder0);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    nodeIndsOverL = find(xyphiarray(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    if any(nodeIndsOverL)
                        xyphiarray(nodeIndsOverL)  = 2*L(dimCtr) - xyphiarray(nodeIndsOverL);
                    end
                else % scalar domain size
                    nodeIndsOverL = find(xyphiarray(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    if any(nodeIndsOverL)
                        xyphiarray(nodeIndsOverL)  = 2*L - xyphiarray(nodeIndsOverL);
                    end
                end
                if any(nodeIndsUnder0)|any(nodeIndsOverL)
                    % change direction of movement upon reflection
                    xyphiarray(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr)) = ... % ugly use of indexing
                        reflectDirection2D(...
                        xyphiarray(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr))...
                        ,dimCtr);
                end
            end
    end
end

end

