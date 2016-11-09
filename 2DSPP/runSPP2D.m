% A function to run 2D SPP model simulations
% LJ Schumacher November 2016
%
% inputs
% T: simulation duration (number of time-steps)
% N: number of objects
% L: size of region containing initial positions
%
% optional inputs
% v0: speed (default 0.05)
% rc: core repulsion radius (default 0.2)
% bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can
%   be single number or 2 element array {'bcx','bcy'} for different
%   bcs along different dimensions
%
% outputs
% xyphiarray: Array containing the position, and movement direction for
% every object and time-point. Format is N by [x,y,phi] by T.

% issues/to-do's:
% - 

function xyphiarray = runSPP2D(T,N,L,varargin)

% parse inputs (see matlab documentation)
iP = inputParser;
checkInt = @(x) x>0&&~mod(x,1);
addRequired(iP,'T',checkInt);
addRequired(iP,'N',checkInt);
addRequired(iP,'L',@checkL);
addOptional(iP,'v0',0.05,@isnumeric)
addOptional(iP,'rc',0.2,@isnumeric)
addOptional(iP,'bc','free',@checkBcs)
parse(iP,T,N,L,varargin{:})
v0 = iP.Results.v0;
rc = iP.Results.rc;
bc = iP.Results.bc;

xyphiarray = NaN(N,3,T);
xyphiarray = initialiseObjects2D(xyphiarray,L);

for t=2:T
    % update direction
    xyphiarray(:,:,t) = updateDirection2D(xyphiarray(:,:,t),xyphiarray(:,:,t-1),L,rc,bc);
    % update position
    xyphiarray(:,:,t) = updatePosition2D(xyphiarray(:,:,t),xyphiarray(:,:,t-1),v0,bc,L);
end
end

function LCheck = checkL(x)
LCheck = false;
if numel(x)==2||numel(x)==1
    LCheck = isnumeric(x);
end
end

function BcCheck = checkBcs(x)
validBcs = {'free','noflux','periodic'};
if iscell(x)&&numel(x)==2
    BcCheck = any(validatestring(x{1},validBcs))...
        &any(validatestring(x{2},validBcs))...
        &any(validatestring(x{3},validBcs));
else
    BcCheck = any(validatestring(x,validBcs));
end
end