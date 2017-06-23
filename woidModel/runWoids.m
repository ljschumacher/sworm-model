% A function to run 2D woid model simulations
% LJ Schumacher November 2016
%
% INPUTS
% T: simulation duration (in seconds)
% N: number of objects
% M: number of nodes in each object
% L: size of region containing initial positions - scalar will give circle
% of radius L, [Lx Ly] will give rectangular domain
%
% optional inputs
% -- general parameters--
% v0: speed (default 0.05)
% dT: time step, scaled adaptively when forces are large (default 1/9s)
% rc: core repulsion radius (default 0.07/2 mm)
% segmentLength: length of a segment between nodes (default 1.2mm/(M-1))
% bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can
%   be single number or 2 element array {'bcx','bcy'} for different
%   bcs along different dimensions
% k_l: stiffness of linear springs connecting nodes
% k_theta: stiffness of rotational springs at nodes
% -- undulation parameters --
% omega_m: angular frequency of undulations (default 0.6Hz)
% theta_0: amplitude of undulations (default pi/4)
% deltaPhase: phase shift in undulations from node to node
% -- reversal parameters --
% revRate: rate for poisson-distributed reversals (default 1/13s)
% revRateCluster: rate for reversals when in a cluster (default 1/130s)
% revRateClusterEdge: increased reversal rates, when worms are poking out of a cluster with their head or tail (default 1/13s no increase)
% revTime: duration of reversal events (default 2s, rounded to integer number of time-steps)
% headNodes: which nodes count as head for defining cluster status, default front 10%
% tailNodes: which nodes count as tail for defining cluster status, default back 10%
% ri: radius at which worms register contact (default 3/2 rc)
% -- slow-down parameters --
% vs: speed when slowed down (default v0/3)
% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
% -- Lennard-Jones parameters --
% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
% eps_LJ: strength of LJ-potential
%
% OUTPUTS
% xyarray: Array containing the position, and movement direction for
% every object and time-point. Format is N by M by [x,y] by T.

% issues/to-do's:
% - is there a more efficient way of storing the coords than 4-d array?
% - make object-orientated, see eg ../woid.m class, vector of woid
% objects...
% - is it still necessary to keep track of node orientation?

function xyarray = runWoids(T,N,M,L,varargin)

% parse inputs (see matlab documentation)
iP = inputParser;
checkInt = @(x) all(x>0&~mod(x,1));
addRequired(iP,'T',checkInt);
addRequired(iP,'N',checkInt);
addRequired(iP,'M',checkInt);
addRequired(iP,'L',@checkL);
addOptional(iP,'v0',0.33,@isnumeric) % worm forward speed is approx 330mu/s
addOptional(iP,'dT',1/9,@isnumeric) % timestep, default 1/9 seconds
addOptional(iP,'rc',0.035,@isnumeric) % worm width is approx 50 to 90 mu = approx 0.07mm
addOptional(iP,'segmentLength',1.2/(M - 1),@isnumeric) % worm length is approx 1.2 mm
addOptional(iP,'bc','free',@checkBcs)
addOptional(iP,'k_l',40,@isnumeric) % stiffness of linear springs connecting nodes
addOptional(iP,'k_theta',20,@isnumeric) % stiffness of rotational springs at nodes
% undulations
addOptional(iP,'omega_m',2*pi*0.6,@isnumeric) % angular frequency of oscillation of movement direction, default 0.6 Hz
addOptional(iP,'theta_0',pi/4,@isnumeric) % amplitude of oscillation of movement direction, default pi/4
addOptional(iP,'deltaPhase',0.24,@isnumeric) % for phase shift in undulations and initial positions, default 0.11
% reversals
addOptional(iP,'revRate',1/13,@isnumeric) % rate for poisson-distributed reversals, default 1/13s
addOptional(iP,'revRateCluster',1/13,@isnumeric) % reduced reversal rates, when worms are in cluster
addOptional(iP,'revRateClusterEdge',1/13,@isnumeric) % increased reversal rates, when worms are poking out of a cluster with their head or tail, default 1/13s (no increase)
addOptional(iP,'revTime',2,@isnumeric) % duration of reversal events, default 2s (will be rounded to integer number of time-steps)
addOptional(iP,'headNodes',1:max(round(M/10),1),checkInt) % which nodes count as head, default front 10%
addOptional(iP,'tailNodes',(M-max(round(M/10),1)+1):M,checkInt) % which nodes count as tail, default back 10%
addOptional(iP,'ri',0.035*3,@isnumeric) % radius at which worms register contact, default 3 rc
% slowing down
addOptional(iP,'vs',0.33/3,@isnumeric) % speed when slowed down, default v0/3
addOptional(iP,'slowingNodes',[1:max(round(M/10),1) (M-max(round(M/10),1)+1):M],checkInt) % which nodes sense proximity, default head and tail
% Lennard-Jones
addOptional(iP,'r_LJcutoff',0,@isnumeric) % cut-off above which lennard jones potential is not acting anymore
addOptional(iP,'eps_LJ',1e-6,@isnumeric) % strength of LJ-potential
addOptional(iP,'sigma_LJ',0,@isnumeric) % particle size for Lennard-Jones force

parse(iP,T,N,M,L,varargin{:})
dT0 = iP.Results.dT;
dT = dT0; % set initial time-step (will be adapted during simulation)
displayOutputEvery = round(1/dT0);
numTimepoints = floor(T/dT0);
v0 = iP.Results.v0;
dTmin = dT0/10/N*v0; % set a mininum timestep below which dT won't be adapted
rc = iP.Results.rc;
bc = iP.Results.bc;
if M>1
    segmentLength = iP.Results.segmentLength;
else
    segmentLength = 0;
end
k_l = iP.Results.k_l;
k_theta = iP.Results.k_theta; % scale with dT as F~v~dT
deltaPhase = iP.Results.deltaPhase;
omega_m = iP.Results.omega_m;
theta_0 = iP.Results.theta_0;
revRate = iP.Results.revRate;
revRateCluster = iP.Results.revRateCluster;
revRateClusterEdge = iP.Results.revRateClusterEdge;
revTime = round(iP.Results.revTime/dT0); % convert to unit of timesteps
headNodes = iP.Results.headNodes;
tailNodes = iP.Results.tailNodes;
ri = iP.Results.ri;
vs = iP.Results.vs;
slowingNodes = iP.Results.slowingNodes;
r_LJcutoff = iP.Results.r_LJcutoff;
eps_LJ = iP.Results.eps_LJ;
sigma_LJ = iP.Results.sigma_LJ;

% check input relationships to each other
% assert(segmentLength>2*rc,...
%     'Segment length must be bigger than node diameter (2*rc). Decrease segment number (M), rc, or increase segmentLength')
assert(min(L)>segmentLength*(M - 1),...
    'Domain size (L) must be bigger than object length (segmentLength*M). Increase L.')
assert(v0>=vs,'vs should be chosen smaller or equal to v0')

% preallocate internal oscillators
theta = NaN(N,M,numTimepoints);
% preallocate reversal states
reversalLogInd = false(N,numTimepoints);
reversalLogIndPrev = reversalLogInd(:,1);
% random phase offset for each object plus phase shift for each node
phaseOffset = wrapTo2Pi(rand(N,1)*2*pi - deltaPhase*(1:M));
% initialise worm positions and node directions - respecting volume
% exclusion
initialExclusionRadius = max(rc,sigma_LJ/2); % so that we don't get too overlapping initial positions, even when rc = 0
[xyarray, theta(:,:,1)] = initialiseWoids(N,M,numTimepoints,L,segmentLength,...
    phaseOffset,theta_0,initialExclusionRadius,bc);
positions = xyarray(:,:,:,1);
orientations = theta(:,:,1);
% initialise time
t = 0;
timeCtr = 1;
disp('Running simulation...')
while t<T
    % find distances between all pairs of objects
    if N==40&&M==49&&numel(L)==2&&~iscell(bc) % check if we can use compiled mex function
        distanceMatrixXY = computeWoidDistancesWithBCs_mex(positions,L,bc);
    elseif N==60&&M==49&&numel(L)==2&&~iscell(bc) % check if we can use compiled mex function
        distanceMatrixXY = computeWoidDistancesWithBCs_N60_mex(positions,L,bc);
    elseif N==80&&M==49&&numel(L)==2&&~iscell(bc) % check if we can use compiled mex function
        distanceMatrixXY = computeWoidDistancesWithBCs_N80_mex(positions,L,bc);
    else
        distanceMatrixXY = computeWoidDistancesWithBCs(positions,L,bc);
    end
    distanceMatrix = sqrt(sum(distanceMatrixXY.^2,5)); % reduce to scalar distances
    % check if any woids are slowed down by neighbors
    [ v, omega ] = slowWorms(distanceMatrix,ri,slowingNodes,vs,v0,omega_m);
    % check if any worms are reversing due to contacts
    reversalLogIndPrev(reversalLogInd(:,timeCtr)) = true; % only update events that happen between timeCtr updates, ie reversal starts
    reversalLogInd = generateReversals(reversalLogInd,timeCtr,distanceMatrix,...
        ri,headNodes,tailNodes,dT,revRate,revTime,revRateCluster,revRateClusterEdge);
    % update internal oscillators / headings
    [orientations, phaseOffset] = updateWoidOscillators(orientations,theta_0,...
        omega,dT,phaseOffset,deltaPhase,[reversalLogIndPrev, reversalLogInd(:, timeCtr)]);
    % calculate forces
    forceArray = calculateForces(distanceMatrixXY,distanceMatrix,...
        rc,orientations,reversalLogInd(:,timeCtr),segmentLength,...
        v,k_l,k_theta*v./v0,theta_0,phaseOffset,sigma_LJ,r_LJcutoff, eps_LJ);
%     % uncomment for debugging...
%     plot(squeeze(positions(objCtr,:,1)),squeeze(positions(objCtr,:,2)),'.-'), axis equal, hold on
%     quiver(squeeze(positions(objCtr,:,1)),squeeze(positions(objCtr,:,2)),...
%         squeeze(forceArray(objCtr,:,1)),squeeze(forceArray(objCtr,:,2)),1)
%     1;
    assert(~any(isinf(forceArray(:))|isnan(forceArray(:))),'Can an unstoppable force move an immovable object? Er...')
    % adapt time-step such that it scales inversily with the max force
    dT = adaptTimeStep(dT0,v0,forceArray);
    if dT<=dTmin
        warning(['Minimum time-step of ' num2str(dTmin) ' reached at time ' num2str(t)])
        dT = dTmin;
    end
    % update position (with boundary conditions)
    [positions, orientations] = applyForces(positions,forceArray,...
        dT,orientations,bc,L);
    assert(~any(isinf(positions(:))),'Uh-oh, something has gone wrong... (infinite positions)')
    % update time
    t = t + dT;
    % output positions and orientations
    if t>=timeCtr*dT0
        reversalLogIndPrev = reversalLogInd(:,timeCtr); % keep this so that we detect end of (fixed-duration) reversals
        timeCtr = timeCtr + 1;
        if mod(timeCtr,displayOutputEvery)==0
            disp(['time = ' num2str(t) ' out of ' num2str(T) 'at ' datestr(now)])
        end
        xyarray(:,:,:,timeCtr) = positions;
        theta(:,:,timeCtr) = orientations;
    end
end
if any(isnan(xyarray(:)))
    warning('some positions remained NaN. this could happen if the adaptive timestep was allowed to increase beyond the default. saving only non-NaNs')
    xyarray = xyarray(~isnan(xyarray(:)));
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