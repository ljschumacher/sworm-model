% A function to run 2D Self-propelled chain model simulations
% LJ Schumacher November 2016
%
% inputs
% T: simulation duration (number of time-steps)
% N: number of objects
% M: number of nodes in each object
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
% every object and time-point. Format is N by M by [x,y,phi] by T.

% issues/to-do's:
% - is there a more efficient way of storing the coords than 4-d array?
% - make object-orientated, see eg ../woid.m class

function xyphiarray = run2DSPC(T,N,M,L,varargin)

% parse inputs (see matlab documentation)
iP = inputParser;
checkInt = @(x) x>0&&~mod(x,1);
addRequired(iP,'T',checkInt);
addRequired(iP,'N',checkInt);
addRequired(iP,'M',checkInt);
addRequired(iP,'L',@checkL);
addOptional(iP,'v0',0.33,@isnumeric) % worm forward speed is approx 330mu/s
addOptional(iP,'dT',1/9,@isnumeric) % adjusts speed, default 1/9 seconds
addOptional(iP,'rc',0.035,@isnumeric) % worm width is approx 50 to 90 mu = approx 0.07mm
addOptional(iP,'segmentLength',1.2/M,@isnumeric) % worm length is approx 1.2 mm
addOptional(iP,'deltaTheta',pi/M,@isnumeric)
addOptional(iP,'bc','free',@checkBcs)
parse(iP,T,N,M,L,varargin{:})
dT = iP.Results.dT;
v0 = iP.Results.v0*dT;
rc = iP.Results.rc;
bc = iP.Results.bc;
segmentLength = iP.Results.segmentLength;
deltaTheta = iP.Results.deltaTheta;

% check input relationships to each other
assert(segmentLength>2*rc,...
    'Segment length must be bigger than node diameter (2*rc). Decrease segment number (M), rc, or increase segmentLength')
assert(min(L)>segmentLength*M,...
    'Domain size (L) must be bigger than object length (segmentLength*M). Increase L.')


xyphiarray = NaN(N,M,3,T);
xyphiarray = initialiseChains2D(xyphiarray,L,segmentLength,deltaTheta);

for t=2:T
    % update direction
    xyphiarray(:,:,:,t) = updateChainDirection2D(xyphiarray(:,:,:,t),xyphiarray(:,:,:,t-1),L,rc,bc);
    % update position
    xyphiarray(:,:,:,t) = updateChainPosition2D(xyphiarray(:,:,:,t),xyphiarray(:,:,:,t-1),v0,bc,L,segmentLength);
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