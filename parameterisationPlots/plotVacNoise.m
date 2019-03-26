% plot vac of simulated worm tracks with different amount of noise

% test simplified woid model with rods

% issues/todo:
% - always reseed random number generator before each simulation?

clear
close all

addpath('../')
addpath('../visualisation')
addpath('../analysis')
addpath('../mex')
% set pdf export options
exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',12,...
    'LineWidth',1);

% general model parameters for all test - unless set otherwise
M = 18; % M: number of nodes in each object
L = [7.5, 7.5];%[20, 20]; % L: size of region containing initial positions - scalar will give circle of radius L, [Lx Ly] will give rectangular domain
param.v0 = 0.33; % v0: speed (default 0.05)
rc = 0.035;
param.rc = 0; % rc: core repulsion radius (default 0.07 mm)
T = 50; % T: simulation duration (number of time-steps)
param.dT = rc/param.v0/16; % dT: time step, gets adapted in simulation
param.saveEvery = round(1/param.dT);
param.bc = 'periodic'; % bc: boundary condition, 'free', 'periodic', or 'noflux' (default 'free'), can be single number or 2 element array {'bcx','bcy'} for different bcs along different dimensions
% undulations
param.omega_m = 0; % angular frequency of oscillation of movement direction, default 0.6 Hz
param.theta_0 = 0; % amplitude of oscillation of movement direction, default pi/4
param.deltaPhase = 0; % for phase shift in undulations and initial positions, default 0.11
% -- reversal parameters --
param.ri = 3*rc;% ri: radius at which rods register contact (default 3 rc)
% -- slow-down parameters --
param.vs = param.v0;% vs: speed when slowed down (default v0/3)
param.slowingNodes = 1:M;% slowingNodes: which nodes register contact (default [1 M], ie head and tail)
% -- Lennard-Jones parameters --
param.r_LJcutoff = 4*rc;% r_LJcutoff: cut-off above which LJ-force is not acting anymore (default 0)
param.sigma_LJ = 0;
param.eps_LJ = 0;% eps_LJ: strength of LJ-potential
food = [];

% test angle noise 
N = 40;
M = 18;
T = round(200); % for N2, scale by *0.33/0.14
maxLag = round(30); % for N2, scale by *0.33/0.14
figure, hold on
noiseLevels = [0.04, 0.05, 0.08]; % for N2, scale by *sqrt(0.14/0.33)
param.ri = 0; % no interaction
param.v0 = 0.33;
param.vs = 0.33;
param.dT = rc/0.33/8; % dT: time step, gets adapted in simulation
param.saveEvery = round(1/2/param.dT);
for noiseLevel = noiseLevels
    param.angleNoise = noiseLevel;% not much point making this any bigger than 10, because it's angular, and even smaller when worm has no stiffness
    param.bc = 'free';
    param.k_theta = 0;
    L = [3 3];
    rng(1)
    xyarray = runWoids(T,N,M,L,param);
    filename = ['vac_M' num2str(M)...
        'angleNoise' num2str(param.angleNoise) '_ktheta_' num2str(param.k_theta)...
        '_dT' num2str(param.dT) '_npr1'];
    % plot velocity autocorrelation of worm worm trajectory to calibrate
    % persistence
    for nn = 1:N
    velocities = gradient(squeeze(xyarray(nn,2,:,:)))./(param.dT*param.saveEvery); % head orientation
    vac(nn,:) = vectorAutoCorrelation(velocities,round(maxLag./(param.dT*param.saveEvery)));
    end
    vac = mean(vac);
    plot(linspace(0,maxLag,length(vac)),vac)
    % plot average heading change vs time ?
end
ylabel('normalised velocity autocorrelation')
xlabel('lag (s)')
h = refline(0,0.23);
h.LineStyle = '--';
h.Color = [0 0 0];
ylim([0 1])
xlim([0 maxLag])
xticks([0 10 20 30])
box on, grid on
lh = legend(num2str(noiseLevels'),'Location','NorthEast');
lh.Title.String = 'noise \eta';
set(gcf,'PaperUnits','centimeters')
exportfig(gcf,[filename '_vac.eps'],exportOptions);
system(['epstopdf ' filename '_vac.eps']);
system(['rm ' filename '_vac.eps']);