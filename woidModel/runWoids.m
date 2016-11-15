% A function to run 2D woid model simulations
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
% - make object-orientated, see eg ../woid.m class, vector of woid
% objects...

function xyphiarray = runWoids(T,N,M,L,varargin)

% parse inputs (see matlab documentation)
iP = inputParser;
checkInt = @(x) x>0&&~mod(x,1);
addRequired(iP,'T',checkInt);
addRequired(iP,'N',checkInt);
addRequired(iP,'M',checkInt);
addRequired(iP,'L',@checkL);
addOptional(iP,'v0',0.33,@isnumeric) % worm forward speed is approx 330mu/s
addOptional(iP,'dT',1/9,@isnumeric) % adjusts speed and undulataions, default 1/9 seconds
addOptional(iP,'rc',0.035,@isnumeric) % worm width is approx 50 to 90 mu = approx 0.07mm
addOptional(iP,'segmentLength',1.2/M,@isnumeric) % worm length is approx 1.2 mm
% undulations
addOptional(iP,'omega_m',2*pi*0.6,@isnumeric) % angular frequency of oscillation of movement direction, default 0.6 Hz
addOptional(iP,'theta_0',pi/4,@isnumeric) % amplitude of oscillation of movement direction, default pi/4
addOptional(iP,'bc','free',@checkBcs)
addOptional(iP,'deltaPhase',2*pi/M*1.2/0.62,@isnumeric) % for phase shift in undulations and initial positions, default given by primary worm wavelength, 620 mu
% reversals
addOptional(iP,'revRate',1/13,@isnumeric) % rate for poisson-distributed reversals, default 1/13s
addOptional(iP,'revTime',2,@numeric) % duration of reversal events, default 2s (will be rounded to integer number of time-steps)

parse(iP,T,N,M,L,varargin{:})
dT = iP.Results.dT;
v0 = iP.Results.v0*dT;
rc = iP.Results.rc;
bc = iP.Results.bc;
segmentLength = iP.Results.segmentLength;
deltaPhase = iP.Results.deltaPhase;
omega_m = iP.Results.omega_m*dT;
theta_0 = iP.Results.theta_0;
revRate = iP.Results.revRate*dT;
revTime = round(iP.Results.revTime/dT);

% check input relationships to each other
assert(segmentLength>2*rc,...
    'Segment length must be bigger than node diameter (2*rc). Decrease segment number (M), rc, or increase segmentLength')
assert(min(L)>segmentLength*M,...
    'Domain size (L) must be bigger than object length (segmentLength*M). Increase L.')


xyphiarray = NaN(N,M,3,T);
% generate internal oscillators 
theta = NaN(N,M,T);
phaseOffset = rand(N,1)*2*pi*ones(1,M) - ones(N,1)*deltaPhase*(1:M); % for each object with random phase offset plus phase shift for each node
theta(:,:,1) = theta_0*sin(omega_m*ones(N,M) + phaseOffset);
% generate reversal times
reversalLogi = logical(poissrnd(revRate,N,T)); % generate reversal events
[reversalWormInd, reversalTimeInd] = find(reversalLogi);
for revCtr = 1:length(reversalWormInd)
    revStart = reversalTimeInd(revCtr);
    if ~reversalLogi(revStart+1) % ignore reversals during reversals
    revEnd = min(revStart + revTime - 1,T);
    reversalLogi(reversalWormInd(revCtr),revStart:revEnd) = true; % set length of all reversal events
    end
end

% initialise worm positions and node directions
xyphiarray = initialiseWoids(xyphiarray,L,segmentLength,theta(:,:,1));

for t=2:T
    % update internal oscillators
    theta(:,:,t) = updateWoidOscillators(theta(:,:,t-1),theta_0,omega_m,t,phaseOffset,reversalLogi(:,t));
    % update direction
    xyphiarray(:,:,:,t) = updateWoidDirection(xyphiarray(:,:,:,t),...
        xyphiarray(:,:,:,t-1),L,rc,bc,theta(:,:,(t - 1):t),reversalLogi(:,t));
    % update position
    xyphiarray(:,:,:,t) = updateWoidPosition(xyphiarray(:,:,:,t),...
        xyphiarray(:,:,:,t-1),v0,bc,L,segmentLength,reversalLogi(:,t));
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