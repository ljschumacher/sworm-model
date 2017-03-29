function [ ] = drawWoidSchematic( )
close all
% draws a schematic of the agent based model
N = 1;
M = 49;
rc = 0.035;
L = 1.25;
segmentLength = 1.2/48;
deltaPhase = 0.25;
theta_0 = pi/4;
bc = 'free';
angles = linspace(0,2*pi,48)'; % for plotting node size

phaseOffset = wrapTo2Pi(deltaPhase*(1:M));
rng(1)
[xyarray, ~] = initialiseWoids(N,M,1,L,segmentLength,phaseOffset,theta_0,rc,bc);

patch(xyarray(:,:,1,1) + rc*cos(angles),...
            xyarray(:,:,2,1) + rc*sin(angles),...
            0.5*ones(1,3),'EdgeColor','None')
axis equal, hold on
plot(xyarray(:,:,1),xyarray(:,:,2),'-',...
    'Marker','.','Color','k','MarkerSize',12);
plot(xyarray(:,1,1),xyarray(:,1,2),'-',...
    'Marker','.','Color','r','MarkerSize',12);
ax = gca;
ax.Color = 'none';
ax.Visible = 'off';

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'LineWidth',1);

set(gcf,'PaperUnits','centimeters')
filename = ['modelSchematic'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

end

