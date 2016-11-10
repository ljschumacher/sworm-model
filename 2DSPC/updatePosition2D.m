function arrayOut = updatePosition2D(arrayNow,arrayPrev,v0,bc,L)
% update positions based on current directions

% issues/to-do's:
%   - could try to optimise the code for the vector domain size/mixed
%   boundary condtion cases to get rid of loops...
%  - mixed periodic boundary conditions can be quite slow

% short-hand for indexing coordinates
x =     1;
y =     2;
phi =   3;

% update position
arrayNow(:,x) = arrayPrev(:,x) + ...
    v0*cos(arrayNow(:,phi));
arrayNow(:,y) = arrayPrev(:,y) + ...
    v0*sin(arrayNow(:,phi));

% check boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can
%   be single number or 2 element array {'bcx','bcy'} for different
%   bcs along different dimensions

if iscell(bc)&&numel(bc)==2
    for dimCtr = [x y]
        switch bc{dimCtr}
            case 'periodic'
                cellIdcsUnder0 = arrayNow(:,dimCtr)<0;
                if numel(L)==3 % vector domain size [L_x L_y]
                    arrayNow(cellIdcsUnder0,dimCtr)  = arrayNow(cellIdcsUnder0,dimCtr) + L(dimCtr);
                    cellIdcsOverL = arrayNow(:,dimCtr)>=L(dimCtr);
                    arrayNow(cellIdcsOverL,dimCtr)  = arrayNow(cellIdcsOverL,dimCtr) - L(dimCtr);
                else % scalar domain size
                    arrayNow(cellIdcsUnder0,dimCtr)  = arrayNow(cellIdcsUnder0,dimCtr) + L;
                    cellIdcsOverL = arrayNow(:,dimCtr)>=L;
                    arrayNow(cellIdcsOverL,dimCtr)  = arrayNow(cellIdcsOverL,dimCtr) - L;
                end
            case 'noflux'
                cellIdcsUnder0 = arrayNow(:,dimCtr)<0;
                arrayNow(cellIdcsUnder0,dimCtr)  = - arrayNow(cellIdcsUnder0,dimCtr);
                if numel(L)==3 % vector domain size [L_x L_y]
                    cellIdcsOverL = arrayNow(:,dimCtr)>=L(dimCtr);
                    arrayNow(cellIdcsOverL,dimCtr)  = 2*L(dimCtr) - arrayNow(cellIdcsOverL,dimCtr);
                else % scalar domain size
                    cellIdcsOverL = arrayNow(:,dimCtr)>=L;
                    arrayNow(cellIdcsOverL,dimCtr)  = 2*L - arrayNow(cellIdcsOverL,dimCtr);
                end
                % change direction of movement upon reflection
                arrayNow(cellIdcsUnder0|cellIdcsOverL,phi) = ...
                    reflectDirection2D(arrayNow(cellIdcsUnder0|cellIdcsOverL,phi),dimCtr);
        end
    end
else
    switch bc
        case 'periodic'
            cellIdcsUnder0 = arrayNow(:,[x y])<0;
            if numel(L)==3 % vector domain size [L_x L_y]
                for dimCtr = [x y]
                    arrayNow(cellIdcsUnder0)  = arrayNow(cellIdcsUnder0) + L(dimCtr);
                    cellIdcsOverL = arrayNow(:,dimCtr)>=L(dimCtr);
                    arrayNow(cellIdcsOverL,dimCtr)  = arrayNow(cellIdcsOverL,dimCtr) - L(dimCtr);
                end
            else % scalar domain size
                arrayNow(cellIdcsUnder0)  = arrayNow(cellIdcsUnder0) + L;
                cellIdcsOverL = arrayNow(:,[x y])>=L;
                arrayNow(cellIdcsOverL)  = arrayNow(cellIdcsOverL) - L;
            end
        case 'noflux'
            for dimCtr = [x y]
                cellIdcsUnder0 = arrayNow(:,dimCtr)<0;
                arrayNow(cellIdcsUnder0,dimCtr)  = - arrayNow(cellIdcsUnder0,dimCtr);
                if numel(L)==3 % vector domain size [L_x L_y]
                    cellIdcsOverL = arrayNow(:,dimCtr)>=L(dimCtr);
                    arrayNow(cellIdcsOverL,dimCtr)  = 2*L(dimCtr) - arrayNow(cellIdcsOverL,dimCtr);
                else % scalar domain size
                    cellIdcsOverL = arrayNow(:,dimCtr)>=L;
                    arrayNow(cellIdcsOverL,dimCtr)  = 2*L - arrayNow(cellIdcsOverL,dimCtr);
                end
                % change direction of movement upon reflection
                arrayNow(cellIdcsUnder0|cellIdcsOverL,phi) = ...
                    reflectDirection2D(arrayNow(cellIdcsUnder0|cellIdcsOverL,phi),dimCtr);
            end
    end
end


arrayOut = arrayNow;