function ax = plotWoidTrajectoriesSingleFrame(xyarray,L,rc,plotColors,centering)
% takes in an array of N by M by x y  by T and makes a movie of the resulting
% trajectories

% issues/to-do:

% short-hand for indexing coordinates
x =     1;
y =     2;

if nargin < 4
    centering = false;
    if nargin <3
        rc = 0.035;
        disp(['Setting node radius to ' num2str(rc) ' for animations.'])
    end
end

N = size(xyarray,1);
M = size(xyarray,2);
if nargin < 4
    plotColors = lines(N);
elseif size(plotColors,1)==1
    plotColors = repmat(plotColors,N,1);
end
assert(size(plotColors,1)==N,'Number of colors not matching number of objects')
angles = linspace(0,2*pi,10)'; % for plotting node size

if centering&&numel(L)==2 % center plot on center of mass - useful for periodic boundary conditions
    xoffset = mean(mean(xyarray(:,:,x),2),1) - L(x)/2;
    yoffset = mean(mean(xyarray(:,:,y),2),1) - L(y)/2;
    xyarray(:,:,x) = xyarray(:,:,x) - xoffset;
    xyarray(:,:,y) = xyarray(:,:,y) - yoffset;
    % re-enforce periodic boundaries
    [ xyarray, ~ ] = checkWoidBoundaryConditions(xyarray, [], 'periodic', L);
end
% set overall axes limits
xrange = minmax(reshape(xyarray(:,:,x),1,numel(xyarray(:,:,x))));
yrange = minmax(reshape(xyarray(:,:,y),1,numel(xyarray(:,:,y))));

% xrange = [floor(xrange(1)) ceil(xrange(2))];
% yrange = [floor(yrange(1)) ceil(yrange(2))];

% if M>1% plot connecting lines btw nodes
%         % don't plot connecting lines for objects that span across a
%         % periodic boundary
%         excludedObjects = any(any(abs(diff(xyarray(:,:,:),1,2))>min(L./2),3),2);
%         plot(xyarray(~excludedObjects,:,x)',xyarray(~excludedObjects,:,y)','k-');
%         if any(~excludedObjects), hold on, end
% elseif rc==0
%         plot(xyarray(:,:,x)',xyarray(:,:,y)','.','Color','k');
%         hold on
% end

ax = gca;
% plot circular domain boundaries, if given scalar domain size
if nargin>=3&&numel(L)==1
    viscircles(ax,[0 0],L,'Color',[0.5 0.5 0.5],'LineWidth',2,'EnhanceVisibility',false);
    hold on
    ax.XLim = [-L L];
    ax.YLim = [-L L];
elseif numel(L)==2
    ax.XLim(1) = min(xrange(1),0);
    ax.XLim(2) = max(xrange(2),L(1));
    ax.YLim(1) = min(yrange(1),0);
    ax.YLim(2) = max(yrange(2),L(2));
else
    ax.XLim = xrange;
    ax.YLim = yrange;
end
ax.DataAspectRatio = [1 1 1];
ax.XTick = [];
ax.YTick = [];
for objCtr = 1:N
    patch(xyarray(objCtr,:,x) + rc*cos(angles),...
        xyarray(objCtr,:,y) + rc*sin(angles),...
        plotColors(objCtr,:),'EdgeColor',plotColors(objCtr,:))
    if objCtr==1
        hold on
    end
    % patches seems to be faster than viscircles
    %         viscircles(squeeze(xyarray(objCtr,:,[x y],frameCtr)),rc*ones(M,1),...
    %             'Color',plotColors(objCtr,:),'EnhanceVisibility',false,...
    %             'LineWidth',1);
end

end