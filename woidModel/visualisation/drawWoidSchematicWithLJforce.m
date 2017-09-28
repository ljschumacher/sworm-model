function [ ] = drawWoidSchematicWithLJforce( )
close all
addpath('../')
% draws a schematic of the agent based model
N = 1;
M = 49;
rc = 0.035;
L = 1.25;
segmentLength = 1.13/(M-1);
deltaPhase = 0.25;
theta_0 = pi/4;
bc = 'free';
angles = linspace(0,2*pi,48)'; % for plotting node size

phaseOffset = wrapTo2Pi(deltaPhase*(1:M));
rng(1)
[xyarray, ~] = initialiseWoids(N,M,1,L,segmentLength,phaseOffset,theta_0,rc,bc);

% plot lennard jones force
eps_LJ = 1;
r_LJcutoff = 5*rc;
sigma_LJ = 2*rc;
f_LJ = @(x) 48*eps_LJ./x.*((sigma_LJ./x).^12 - 1/2*(sigma_LJ./x).^6);

x = linspace(-0.1,1.2,2600);
y = linspace(-0.2,0.3,1000);
[X,Y] = meshgrid(x,y);
F_LJ = zeros(size(X));
for nodeCtr = 1:M
    r = sqrt((x - xyarray(:,nodeCtr,1)).^2 + (y - xyarray(:,nodeCtr,2))'.^2);
    % ignore values smaller than the node radius for plotting
    r(r<2*rc) = 2*rc;
    r(r>r_LJcutoff) = Inf;
    F_LJ = F_LJ + f_LJ(r);
end
F_LJ(F_LJ>0) = -min(min(F_LJ));
h = pcolor(X,Y,F_LJ);
h.MeshStyle = 'none';
hold on
% load divergent colormap
divcmap = load('/Users/linus/Dropbox/Utilities/colormaps_ascii/diverging/cmap_BuRd.txt');
colormap(divcmap)

% plot worm
patch(xyarray(:,:,1,1) + rc*cos(angles),...
    xyarray(:,:,2,1) + rc*sin(angles),...
    0.5*ones(1,3),'EdgeColor','None')
plot(xyarray(:,:,1),xyarray(:,:,2),'-',...
    'Marker','.','Color','k','MarkerSize',12);
plot(xyarray(:,1,1),xyarray(:,1,2),'-',...
    'Marker','.','Color','r','MarkerSize',12);
ax = gca;
ax.Color = 'none';
ax.DataAspectRatio = [1 1 1];
ax.Visible = 'on';
ax.XLabel.String = 'x';
ax.YLabel.String = 'y';

hcb = colorbar;
hcb.Label.String = 'Lennard-Jones Force';

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',14,...
    'Resolution',300,...
    'LineWidth',1);

set(gcf,'PaperUnits','centimeters')
filename = ['modelSchematicWithLJForce'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

end

