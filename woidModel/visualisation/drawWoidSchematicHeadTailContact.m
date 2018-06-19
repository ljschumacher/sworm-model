function [ ] = drawWoidSchematicHeadTailContact( )
close all
addpath('../')
% draws a schematic of the agent based model
N = 1;
M = 36;
rc = 0.035;
L = 1.25;
segmentLength = 1.2/(M-1);
deltaPhase = 0.25*49/M;
theta_0 = pi/4;
bc = 'free';
angles = linspace(0,2*pi,48)'; % for plotting node size
ri = 3*rc;

headNodes = 1:max(round(M/10),1);
tailNodes = (M-max(round(M/10),1)+1):M;

phaseOffset = wrapTo2Pi(deltaPhase*(1:M));
rng(1)
[xyarray, ~] = initialiseWoids(N,M,1,L,segmentLength,phaseOffset,theta_0,rc,bc);
%% head green, tail red
figure
% draw head zone
patch(xyarray(:,headNodes,1,1) + ri*cos(angles),...
    xyarray(:,headNodes,2,1) + ri*sin(angles),...
    [0.2 0.8 0.2],'EdgeColor','None')
axis equal, hold on
% draw tail zone
patch(xyarray(:,tailNodes,1,1) + ri*cos(angles),...
    xyarray(:,tailNodes,2,1) + ri*sin(angles),...
    [0.8 0.2 0.2],'EdgeColor','None')
% draw worm
patch(xyarray(:,:,1,1) + rc*cos(angles),...
    xyarray(:,:,2,1) + rc*sin(angles),...
    0.5*ones(1,3),'EdgeColor','None')

ax = gca;
ax.Color = 'none';
ax.Visible = 'off';

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'LineWidth',1);

set(gcf,'PaperUnits','centimeters')
filename = ['modelSchematicHeadContact'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

%% head red, tail green
figure
% draw head zone
patch(xyarray(:,headNodes,1,1) + ri*cos(angles),...
    xyarray(:,headNodes,2,1) + ri*sin(angles),...
    [0.8 0.2 0.2],'EdgeColor','None')
axis equal, hold on
% draw tail zone
patch(xyarray(:,tailNodes,1,1) + ri*cos(angles),...
    xyarray(:,tailNodes,2,1) + ri*sin(angles),...
    [0.2 0.8 0.2],'EdgeColor','None')
% draw worm
patch(xyarray(:,:,1,1) + rc*cos(angles),...
    xyarray(:,:,2,1) + rc*sin(angles),...
    0.5*ones(1,3),'EdgeColor','None')

ax = gca;
ax.Color = 'none';
ax.Visible = 'off';

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'LineWidth',1);

set(gcf,'PaperUnits','centimeters')
filename = ['modelSchematicTailContact'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

%% head green, tail green
figure
% draw head zone
patch(xyarray(:,headNodes,1,1) + ri*cos(angles),...
    xyarray(:,headNodes,2,1) + ri*sin(angles),...
    [0.2 0.8 0.2],'EdgeColor','None')
axis equal, hold on
% draw tail zone
patch(xyarray(:,tailNodes,1,1) + ri*cos(angles),...
    xyarray(:,tailNodes,2,1) + ri*sin(angles),...
    [0.2 0.8 0.2],'EdgeColor','None')
% draw worm
patch(xyarray(:,:,1,1) + rc*cos(angles),...
    xyarray(:,:,2,1) + rc*sin(angles),...
    0.5*ones(1,3),'EdgeColor','None')

ax = gca;
ax.Color = 'none';
ax.Visible = 'off';

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'LineWidth',1);

set(gcf,'PaperUnits','centimeters')
filename = ['modelSchematicBothContact'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

%% head red, tail red
figure
% draw head zone
patch(xyarray(:,headNodes,1,1) + ri*cos(angles),...
    xyarray(:,headNodes,2,1) + ri*sin(angles),...
    [0.8 0.2 0.2],'EdgeColor','None')
axis equal, hold on
% draw tail zone
patch(xyarray(:,tailNodes,1,1) + ri*cos(angles),...
    xyarray(:,tailNodes,2,1) + ri*sin(angles),...
    [0.8 0.2 0.2],'EdgeColor','None')
% draw worm
patch(xyarray(:,:,1,1) + rc*cos(angles),...
    xyarray(:,:,2,1) + rc*sin(angles),...
    0.5*ones(1,3),'EdgeColor','None')

ax = gca;
ax.Color = 'none';
ax.Visible = 'off';

exportOptions = struct('Format','eps2',...
    'Color','rgb',...
    'Width',10,...
    'Resolution',300,...
    'LineWidth',1);

set(gcf,'PaperUnits','centimeters')
filename = ['modelSchematicNoContact'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['rm ' filename '.eps']);

end

