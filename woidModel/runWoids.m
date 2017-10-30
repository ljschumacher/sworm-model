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
% v0: speed (default 0.33 for npr1, 0.14 for N2)
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
% angleNoise: noise (gaussian) in direction of movement for head node
% -- reversal parameters --
% revRate: rate for poisson-distributed reversals (default 1/13s)
% revRateCluster: rate for reversals when in a cluster (default 1/130s)
% revRateClusterEdge: increased reversal rates, when worms are poking out of a cluster with their head or tail (default 1/13s no increase)
% revTime: duration of reversal events (default 2s, rounded to integer number of time-steps)
% headNodes: which nodes count as head for defining cluster status, default front 10%
% tailNodes: which nodes count as tail for defining cluster status, default back 10%
% ri: radius at which worms register contact (default 3/2 rc)
% Rir: relative interaction radius for reversals, default 1 (ie = ri)
% num_nbr_max_per_node: how many max nbrs per node to count for density-dependent speed
% -- slow-down parameters --
% vs: speed when slowed down (default v0/3)
%   for stochastic slowing, low dwelling speeds are 0.018 and 0.014 for
%   npr1/N2
% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
% slowingMode: how worms slow down, 'abrupt', 'gradual', 'density', or
% 'stochastic'
% k_dwell: rate to enter low-speed dwelling state, for stochastic slowing
%   mode, and in the absence of neighbours (default 1/275s for npr1,
%   1/4 s for N2) ???
% dkdN_dwell: increase in dwelling rate with neighbour density
% k_undwell: rate to leave low-speed dwelling state, for stochastic slowing
%   mode, and in the absence of neighbours (default 1/0.9s for npr1,
%   1/2.2s for N2)
% Ris: relative interaction radius for slowing, default 1 (ie = ri)
% -- adhesion (Lennard-Jones) parameters --
% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
% eps_LJ: strength of LJ-potential
% sigma_LJ: particle size of the LJ-force, so minimum is at ca 1.122sigma
% LJnodes: which nodes the LJ-forces acts on (default 1:M). This is
%   asymmetric in the sense that these nodes will feel the force from all
%   other nodes, but the other nodes won't feel the force from these nodes
% -- state parameters --
% k_roam: rate to spontaneously enter the roaming state (like off-food)
% k_unroam: rate to spontaneously leave the roaming state
%
% OUTPUTS
% xyarray: Array containing the position, and movement direction for
% every object and time-point. Format is N by M by [x,y] by T.

% issues/to-do's:
% - is there a more efficient way of storing the coords than 4-d array?
% - make object-orientated, see eg ../woid.m class, vector of woid
% objects...

function [xyarray, currentState] = runWoids(T,N,M,L,varargin)

% parse inputs (see matlab documentation)
iP = inputParser;
checkInt = @(x) all(x>0&~mod(x,1));
isPositiveNumeric = @(x) all(isnumeric(x)&x>0);
addRequired(iP,'T',checkInt);
addRequired(iP,'N',checkInt);
addRequired(iP,'M',checkInt);
addRequired(iP,'L',@checkL);
addOptional(iP,'v0',0.33,@isnumeric) % worm forward speed is approx 330mu/s
addOptional(iP,'dT',0.035/0.33/16,isPositiveNumeric) % timestep, default rc/v0/16 seconds
addOptional(iP,'rc',0.035,@isnumeric) % worm width is approx 50 to 90 mu = approx 0.07mm
addOptional(iP,'segmentLength',(1.2-2*0.035)/(M - 1),@isnumeric) % worm length is approx 1.2 mm
addOptional(iP,'bc','free',@checkBcs)
addOptional(iP,'k_l',40,@isnumeric) % stiffness of linear springs connecting nodes
addOptional(iP,'k_theta',20,@isnumeric) % stiffness of rotational springs at nodes
% undulations
addOptional(iP,'omega_m',2*pi*0.6,@isnumeric) % angular frequency of oscillation of movement direction, default 0.6 Hz
addOptional(iP,'theta_0',pi/4,@isnumeric) % amplitude of oscillation of movement direction, default pi/4
addOptional(iP,'deltaPhase',0.24,@isnumeric) % for phase shift in undulations and initial positions, default 0.11
addOptional(iP,'angleNoise',0,@isnumeric) % noise in heading of worm
% reversals
addOptional(iP,'revRate',0,@isnumeric) % rate for poisson-distributed reversals, default 1/13s
addOptional(iP,'revRateCluster',0,@isnumeric) % reduced reversal rates, when worms are in cluster
addOptional(iP,'revRateClusterEdge',1/13,@isnumeric) % increased reversal rates, when worms are poking out of a cluster with their head or tail, default 1/13s (no increase)
addOptional(iP,'revTime',2,@isnumeric) % duration of reversal events, default 2s (will be rounded to integer number of time-steps)
addOptional(iP,'headNodes',1:max(round(M/10),1),checkInt) % which nodes count as head, default front 10%
addOptional(iP,'tailNodes',(M-max(round(M/10),1)+1):M,checkInt) % which nodes count as tail, default back 10%
addOptional(iP,'ri',0.035*3,@isnumeric) % radius at which worms register contact, default 3 rc
addOptional(iP,'Rir',1,@isnumeric) % relative interaction radius for reversals, default 1 (ie = ri)
% slowing down
addOptional(iP,'vs',0.33/3,@isnumeric) % speed when slowed down, default v0/3
addOptional(iP,'slowingNodes',1:M,checkInt) % which nodes sense proximity, default all
addOptional(iP,'slowingMode','gradual',@checkSlowingMode) % 'gradual' or 'abrupt'or 'density'
addOptional(iP,'Ris',1,@isnumeric) % relative interaction radius for slowing, default 1 (ie = ri)
addOptional(iP,'num_nbr_max_per_node',1,checkInt) % how many max nbrs per node to count for density-dependent speed
addOptional(iP,'k_dwell',0,@isnumeric) % rate to enter low-speed dwelling state, for stochastic slowing
addOptional(iP,'k_undwell',0,@isnumeric) % rate to leave low-speed dwelling state, for stochastic slowing
addOptional(iP,'dkdN_dwell',0,@isnumeric) % increase in dwelling rate with nbr density, for stochastic slowing
% adhesion forces (Lennard-Jones)
addOptional(iP,'r_LJcutoff',0,@isnumeric) % cut-off above which lennard jones potential is not acting anymore
addOptional(iP,'eps_LJ',1e-6,@isnumeric) % strength of LJ-potential
addOptional(iP,'sigma_LJ',0,@isnumeric) % particle size for Lennard-Jones force
addOptional(iP,'LJnodes',1:M,checkInt) % nodes which feel LJ-force
% state parameters
addOptional(iP,'k_roam',0,@isnumeric) % rate to spontaneously enter the roaming state (like off-food) (default 0)
addOptional(iP,'k_unroam',0,@isnumeric) % rate to spontaneously leave the roaming state (default 0)
% simulation parameters
addOptional(iP,'saveEvery',1,checkInt);
addOptional(iP,'resumeState',struct([]),@checkResumeState); % current state of previous simulation to initialize from

parse(iP,T,N,M,L,varargin{:})
dT0 = iP.Results.dT;
dT = dT0; % set initial time-step (will be adapted during simulation)
saveEvery = iP.Results.saveEvery;
displayOutputEvery = round(1/dT0);
numTimepoints = floor(T/dT0);
numSavepoints = floor(T/dT0/saveEvery);
resumeState = iP.Results.resumeState;

v0 = iP.Results.v0;
dTmin = dT0/10/sqrt(N)*v0; % set a mininum timestep below which dT won't be adapted
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
angleNoise = iP.Results.angleNoise;
revRate = iP.Results.revRate;
revRateCluster = iP.Results.revRateCluster;
revRateClusterEdge = iP.Results.revRateClusterEdge;
revTime = round(iP.Results.revTime/dT0); % convert to unit of timesteps
headNodes = iP.Results.headNodes;
tailNodes = iP.Results.tailNodes;
ri = iP.Results.ri;
Rir = iP.Results.Rir;
Ris = iP.Results.Ris;
vs = iP.Results.vs;
slowingNodes = iP.Results.slowingNodes;
slowingMode = iP.Results.slowingMode;
num_nbr_max_per_node = iP.Results.num_nbr_max_per_node;
k_dwell = iP.Results.k_dwell;
k_undwell = iP.Results.k_undwell;
dkdN_dwell = iP.Results.dkdN_dwell;
r_LJcutoff = iP.Results.r_LJcutoff;
eps_LJ = iP.Results.eps_LJ;
sigma_LJ = iP.Results.sigma_LJ;
LJnodes = iP.Results.LJnodes;
k_roam = iP.Results.k_roam;
k_unroam = iP.Results.k_unroam;

% check input relationships to each other
% assert(segmentLength>2*rc,...
%     'Segment length must be bigger than node diameter (2*rc). Decrease segment number (M), rc, or increase segmentLength')
assert(min(L)>segmentLength*(M - 1),...
    'Domain size (L) must be bigger than object length (segmentLength*M). Increase L.')
assert(v0>=vs,'vs should be chosen smaller or equal to v0')

% preallocate reversal states
reversalLogInd = false(N,numTimepoints);
% check if resuming a previous simulation
if isempty(resumeState) 
    % preallocate roaming state variable
    roamingLogInd = false(N,1);
    dwellLogInd = false(N,1);
    % generate random phase offset for each object plus phase shift for each node
    phaseOffset = wrapTo2Pi(rand(N,1)*2*pi - deltaPhase*(1:M));
    % initialise worm positions and node directions - respecting volume exclusion
    initialExclusionRadius = max(rc,sigma_LJ/2); % so that we don't get too overlapping initial positions, even when rc = 0
    [xyarray, orientations] = initialiseWoids(N,M,numSavepoints,L,segmentLength,...
        phaseOffset,theta_0,initialExclusionRadius,bc);
    positions = xyarray(:,:,:,1);
else
    reversalLogInd(:,1) = resumeState.reversalLogInd;
    roamingLogInd = resumeState.roamingLogInd;
    dwellLogInd = resumeState.dwellLogInd;
    phaseOffset = resumeState.phaseOffset;
    positions = resumeState.positions;
    orientations = resumeState.orientations;
    xyarray = NaN(N,M,2,numSavepoints);
    xyarray(:,:,:,1) = positions;
end
reversalLogIndPrev = reversalLogInd(:,1);

% initialise time
t = 0;
timeCtr = 1;
saveCtr = 1;
disp(['Running simulation...' datestr(now)])
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
    % check if any woids are changing their movement state
    roamingLogInd = updateRoamingState(roamingLogInd,k_roam,k_unroam,dT);
    % check if any woids are slowed down by neighbors
    [ v, omega, dwellLogInd ] = slowWorms(distanceMatrix,Ris*ri,slowingNodes,slowingMode,...
        vs,v0,omega_m,num_nbr_max_per_node,roamingLogInd,k_dwell,k_undwell,dkdN_dwell,dwellLogInd,dT);
    % check if any worms are reversing due to contacts
    reversalLogIndPrev(reversalLogInd(:,timeCtr)) = true; % only update events that happen between timeCtr updates, ie reversal starts
    reversalLogInd = generateReversals(reversalLogInd,timeCtr,distanceMatrix,...
        Rir*ri,headNodes,tailNodes,dT,revRate,revTime,revRateCluster,revRateClusterEdge,roamingLogInd);
    % update internal oscillators / headings
    [orientations, phaseOffset] = updateWoidOscillators(orientations,theta_0,...
        omega,dT,phaseOffset,deltaPhase,[reversalLogIndPrev, reversalLogInd(:, timeCtr)]);
    % calculate forces
    forceArray = calculateForces(distanceMatrixXY,distanceMatrix,...
        rc,orientations,reversalLogInd(:,timeCtr),segmentLength,...
        v,k_l,k_theta*v./v0,theta_0,phaseOffset,sigma_LJ,r_LJcutoff,eps_LJ,LJnodes,...
        angleNoise);
    try
    assert(~any(isinf(forceArray(:))|isnan(forceArray(:))),'Can an unstoppable force move an immovable object? Er...')
    catch
        1;
    end
    % adapt time-step such that it scales inversily with the max force
    dT = adaptTimeStep(dT0,v0,forceArray);
    if dT<=dTmin
        warning(['Minimum time-step of ' num2str(dTmin) ' reached at time ' num2str(t)])
        dT = dTmin;
    end
    assert(dT<=dT0)
    % update position (with boundary conditions)
    [positions, orientations] = applyForces(positions,forceArray,...
        dT,orientations,bc,L);
    assert(~any(isinf(positions(:))),'Uh-oh, something has gone wrong... (infinite positions)')
    % update time
    t = t + dT;
    % output positions and orientations
    if t>=timeCtr*dT0
        reversalLogIndPrev = reversalLogInd(:,timeCtr); % keep this so that we detect end of (fixed-duration) reversals
        if t>=saveCtr*dT0*saveEvery
            saveCtr = saveCtr+1;
            xyarray(:,:,:,saveCtr) = positions;
            % save other outputs to enable resuming simulations
            currentState.positions = positions;
            currentState.orientations = orientations;
            currentState.phaseOffset = phaseOffset;
            currentState.reversalLogInd = reversalLogInd(:,timeCtr);
            currentState.roamingLogInd = roamingLogInd;
            currentState.dwellLogInd= dwellLogInd;
        end
        timeCtr = timeCtr + 1;
        if mod(timeCtr,displayOutputEvery)==0
            disp(['time = ' num2str(t) ' out of ' num2str(T) ' at ' datestr(now)])
        end
    end
end
if any(isnan(xyarray(:)))
    warning('some positions remained NaN. this could happen if the adaptive timestep was allowed to increase beyond the default. saving only non-NaNs')
    xyarray = reshape(xyarray(~isnan(xyarray(:))),N,M,2,[]);
end
end

function LCheck = checkL(x)
LCheck = false;
if numel(x)==2||numel(x)==1
    LCheck = isnumeric(x);
end
end

function BcCheck = checkBcs(s)
validBcs = {'free','noflux','periodic'};
if iscell(s)&&numel(s)==2
    BcCheck = any(strcmp(s{1},validBcs))...
        &any(strcmp(s{2},validBcs));
else
    BcCheck = any(validatestring(s,validBcs));
end
end

function SlowModeCheck = checkSlowingMode(s)
validSlowingModes = {'gradual','abrupt','density','stochastic','stochastic_bynode'};
SlowModeCheck = any(strcmp(s,validSlowingModes));
end

function resumeStateCheck = checkResumeState(resumeState)
    resumeStateCheck = isstruct(resumeState)&&...
        all(isfield(resumeState,{'positions','orientations','phaseOffset',...
        'reversalLogInd','roamingLogInd','dwellLogInd'}));
end
