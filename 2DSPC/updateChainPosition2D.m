function arrayOut = updateChainPosition2D(arrayNow,arrayPrev,v0,bc,L)
% update positions based on current directions

% issues/to-do's:
%   - could try to optimise the code for the vector domain size/mixed
%   boundary condition cases to get rid of loops...
%   - mixed periodic boundary conditions can be quite slow?
%   - move check of boundary conditions to separate function?

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

% update position
arrayNow(:,:,x) = arrayPrev(:,:,x) + ...
    v0*cos(arrayNow(:,:,phi));
arrayNow(:,:,y) = arrayPrev(:,:,y) + ...
    v0*sin(arrayNow(:,:,phi));

N = size(arrayNow,1);
M = size(arrayNow,2);
ndim = size(arrayNow,3)-1;

% check boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can
%   be single number or 2 element array {'bcx','bcy'} for different
%   bcs along different dimensions

if iscell(bc)&&numel(bc)==ndim
    for dimCtr = [x y]
        switch bc{dimCtr}
            case 'periodic'
                nodeIndsUnder0 = find(arrayNow(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    arrayNow(nodeIndsUnder0)  = arrayNow(nodeIndsUnder0) + L(dimCtr);
                    nodeIndsOverL = find(arrayNow(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    arrayNow(nodeIndsOverL)  = arrayNow(nodeIndsOverL) - L(dimCtr);
                else % scalar domain size
                    arrayNow(nodeIndsUnder0)  = arrayNow(nodeIndsUnder0) + L;
                    nodeIndsOverL = find(arrayNow(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    arrayNow(nodeIndsOverL)  = arrayNow(nodeIndsOverL) - L;
                end
            case 'noflux'
                nodeIndsUnder0 = find(arrayNow(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                arrayNow(nodeIndsUnder0)  = - arrayNow(nodeIndsUnder0);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    nodeIndsOverL = find(arrayNow(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    arrayNow(nodeIndsOverL)  = 2*L(dimCtr) - arrayNow(nodeIndsOverL);
                else % scalar domain size
                    nodeIndsOverL = find(arrayNow(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    arrayNow(nodeIndsOverL)  = 2*L - arrayNow(nodeIndsOverL);
                end
                % change direction of movement upon reflection
                arrayNow(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr)) = ... % ugly use of indexing
                    reflectDirection2D(...
                    arrayNow(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr))...
                    ,dimCtr);
        end
    end
else
    switch bc
        case 'periodic'
            if numel(L)==ndim % vector domain size [L_x L_y]
                for dimCtr = [x y]
                    nodeIndsUnder0 = find(arrayNow(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                    arrayNow(nodeIndsUnder0)  = arrayNow(nodeIndsUnder0) + L(dimCtr);
                    nodeIndsOverL = find(arrayNow(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    arrayNow(nodeIndsOverL)  = arrayNow(nodeIndsOverL) - L(dimCtr);
                end
            else % scalar domain size
                nodeLogiUnder0 = arrayNow(:,:,[x y])<0;
                arrayNow(nodeLogiUnder0)  = arrayNow(nodeLogiUnder0) + L;
                nodeLogiOverL = arrayNow(:,:,[x y])>=L;
                arrayNow(nodeLogiOverL)  = arrayNow(nodeLogiOverL) - L;
            end
        case 'noflux'
            for dimCtr = [x y]
                nodeIndsUnder0 = find(arrayNow(:,:,dimCtr)<0) + N*M*(dimCtr - 1);
                    arrayNow(nodeIndsUnder0)  = - arrayNow(nodeIndsUnder0);
                if numel(L)==ndim % vector domain size [L_x L_y]
                    nodeIndsOverL = find(arrayNow(:,:,dimCtr)>=L(dimCtr)) + N*M*(dimCtr - 1);
                    if any(nodeIndsOverL)
                        arrayNow(nodeIndsOverL)  = 2*L(dimCtr) - arrayNow(nodeIndsOverL);
                    end
                else % scalar domain size
                    nodeIndsOverL = find(arrayNow(:,:,dimCtr)>=L) + N*M*(dimCtr - 1);
                    if any(nodeIndsOverL)
                        arrayNow(nodeIndsOverL)  = 2*L - arrayNow(nodeIndsOverL);
                    end
                end
                % change direction of movement upon reflection
                arrayNow(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr)) = ... % ugly use of indexing
                    reflectDirection2D(...
                    arrayNow(union(nodeIndsUnder0,nodeIndsOverL) + N*M*(phi - dimCtr))...
                    ,dimCtr);
            end
    end
end


arrayOut = arrayNow;