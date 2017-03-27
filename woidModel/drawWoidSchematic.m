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

phaseOffset = wrapTo2Pi(deltaPhase*(1:M));
[xyarray, theta] = initialiseWoids(N,M,1,L,segmentLength,phaseOffset,theta_0,rc,bc);

viscircles(squeeze(xyarray),rc*ones(M,1),'Color',0.5*ones(3,1),'EnhanceVisibility',false,...
    'LineWidth',1);
axis equal, hold on
plot(xyarray(:,:,1),xyarray(:,:,2),'-',...
    'Marker','.','Color','k');
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

