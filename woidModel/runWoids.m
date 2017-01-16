% A function to run 2D woid model simulations
% LJ Schumacher November 2016
%
% INPUTS
% T: simulation duration (number of time-steps)
% N: number of objects
% M: number of nodes in each object
% L: size of region containing initial positions - scalar will give circle
% of radius L, [Lx Ly] will give rectangular domain
%
% optional inputs
% -- general parameters--
% v0: speed (default 0.05)
% dT: time step, scales other parameters such as velocities and rates
% (default 1/9s)
% rc: core repulsion radius (default 0.07 microns)
% segmentLength: length of a segment between nodes (default 1.2mm/(M-1))
% bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can
%   be single number or 2 element array {'bcx','bcy'} for different
%   bcs along different dimensions
% kl: stiffness of linear springs connecting nodes
% k_theta: stiffness of rotational springs at nodes
% -- undulation parameters --
% omega_m: angular frequency of undulations (default 0.6Hz)
% theta_0: amplitude of undulations (default pi/4)
% deltaPhase: phase shift in undulations from node to node
% -- reversal parameters --
% revRate: rate for poisson-distributed reversals (default 1/13s)
% revTime: duration of reversal events (default 2s, rounded to integer number of time-steps)
% -- slow-down parameters --
% rs: radius at which worms register contact (default 2 rc)
% vs: speed when slowed down (default v0/3)
% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
%
% OUTPUTS
% xyphiarray: Array containing the position, and movement direction for
% every object and time-point. Format is N by M by [x,y,phi] by T.

% issues/to-do's:
% - is there a more efficient way of storing the coords than 4-d array?
% - make object-orientated, see eg ../woid.m class, vector of woid
% objects...
% - is it still necessary to keep track of node orientation?

function xyphiarray = runWoids(T,N,M,L,varargin)

% parse inputs (see matlab documentation)
iP = inputParser;
checkInt = @(x) all(x>0&~mod(x,1));
addRequired(iP,'T',checkInt);
addRequired(iP,'N',checkInt);
addRequired(iP,'M',checkInt);
addRequired(iP,'L',@checkL);
addOptional(iP,'v0',0.33,@isnumeric) % worm forward speed is approx 330mu/s
addOptional(iP,'dT',1/9,@isnumeric) % adjusts speed and undulataions, default 1/9 seconds
addOptional(iP,'rc',0.035,@isnumeric) % worm width is approx 50 to 90 mu = approx 0.07mm
addOptional(iP,'segmentLength',1.2/(M - 1),@isnumeric) % worm length is approx 1.2 mm
addOptional(iP,'bc','free',@checkBcs)
addOptional(iP,'kl',10,@isnumeric) % stiffness of linear springs connecting nodes
addOptional(iP,'k_theta',20,@isnumeric) % stiffness of rotational springs at nodes
% undulations
addOptional(iP,'omega_m',2*pi*0.6,@isnumeric) % angular frequency of oscillation of movement direction, default 0.6 Hz
addOptional(iP,'theta_0',pi/4,@isnumeric) % amplitude of oscillation of movement direction, default pi/4
addOptional(iP,'deltaPhase',0.11,@isnumeric) % for phase shift in undulations and initial positions, default 0.11
% reversals
addOptional(iP,'revRate',1/13,@isnumeric) % rate for poisson-distributed reversals, default 1/13s
addOptional(iP,'revRateCluster',1/130,@isnumeric) % reduced reversal rates, when worms are in cluster
addOptional(iP,'revTime',2,@isnumeric) % duration of reversal events, default 2s (will be rounded to integer number of time-steps)
addOptional(iP,'headNodes',1:max(round(M/10),1),checkInt) % which nodes count as head, default fron 10%
addOptional(iP,'tailNodes',(M-max(round(M/10),1)+1):M,checkInt) % which nodes count as tail, default back 10%
% slowing down
addOptional(iP,'rs',0.035*2,@isnumeric) % radius at which worms slow down, default 2 rc
addOptional(iP,'vs',0.33/3,@isnumeric) % speed when slowed down, default v0/3
addOptional(iP,'slowingNodes',[1:max(round(M/10),1) (M-max(round(M/10),1)+1):M],checkInt) % which nodes sense proximity, default head and tail

parse(iP,T,N,M,L,varargin{:})
dT = iP.Results.dT;
v0 = iP.Results.v0*dT;
rc = iP.Results.rc;
bc = iP.Results.bc;
segmentLength = iP.Results.segmentLength;
kl = iP.Results.kl*dT; % scale with dT as F~v~dT
k_theta = iP.Results.k_theta*dT; % scale with dT as F~v~dT
deltaPhase = iP.Results.deltaPhase;
omega_m = iP.Results.omega_m*dT;
theta_0 = iP.Results.theta_0;
revRate = iP.Results.revRate*dT;
revRateCluster = iP.Results.revRateCluster*dT;
revTime = round(iP.Results.revTime/dT);
headNodes = iP.Results.headNodes;
tailNodes = iP.Results.tailNodes;
rs = iP.Results.rs;
vs = iP.Results.vs*dT;
slowingNodes = iP.Results.slowingNodes;

% check input relationships to each other
% assert(segmentLength>2*rc,...
%     'Segment length must be bigger than node diameter (2*rc). Decrease segment number (M), rc, or increase segmentLength')
assert(min(L)>segmentLength*(M - 1),...
    'Domain size (L) must be bigger than object length (segmentLength*M). Increase L.')
assert(v0>=vs,'vs should be chosen smaller or equal to v0')

xyphiarray = NaN(N,M,3,T);
% generate internal oscillators 
theta = NaN(N,M,T);
phaseOffset = rand(N,1)*2*pi*ones(1,M) - ones(N,1)*deltaPhase*(1:M); % for each object with random phase offset plus phase shift for each node
theta(:,:,1) = theta_0*sin(omega_m*ones(N,M) + phaseOffset);
% preallocate reversal states
reversalLogInd = false(N,T);

% initialise worm positions and node directions - respecting volume
% exclusion
xyphiarray = initialiseWoids(xyphiarray,L,segmentLength,theta(:,:,1),rc,bc);
disp('Running simulation...')
for t=2:T
    % find distances between all pairs of objects
    distanceMatrixXY = computeWoidDistancesWithBCs(xyphiarray(:,:,1:2,t-1),L,bc);
    distanceMatrix = sqrt(sum(distanceMatrixXY.^2,5)); % reduce to scalar
    % check if any woids are slowed down by neighbors
    slowLogInd = findWoidNeighbors(distanceMatrix,2*rs,slowingNodes);
    v = vs*slowLogInd + v0*(~slowLogInd); % adjust speed for slowed worms
    omega = omega_m*(~slowLogInd + vs/v0*slowLogInd); % adjust internal oscillator freq for slowed worms
    % check if any worms are reversing due to contacts
    reversalLogInd = generateReversals(reversalLogInd,t,distanceMatrix,...
    2*rs,headNodes,tailNodes,revRate,revTime,revRateCluster);
    % update internal oscillators
    theta(:,:,t) = updateWoidOscillators(theta(:,:,t-1),theta_0,omega,t,phaseOffset,reversalLogInd(:,t));     
    % calculate forces
    forceArray = calculateForces(xyphiarray(:,:,:,t-1),rc,distanceMatrixXY,...
        distanceMatrix,theta(:,:,(t-1):t),reversalLogInd(:,(t-1):t),segmentLength,v,kl,k_theta);
    assert(~any(isinf(forceArray(:))|isnan(forceArray(:))),'Can an unstoppable force move an immovable object? Er...')
    % update position (with boundary conditions)
    xyphiarray(:,:,:,t) = applyForces(xyphiarray(:,:,:,t-1),forceArray,bc,L);
    assert(~any(isinf(xyphiarray(:))),'Uh-oh, something has gone wrong... (try using a smaller time-step?)')
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