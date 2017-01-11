function [ ] = animateWoidTrajectories(xyphiarray,filename,L,rc)
% takes in an array of N by M by x y phi by T and makes a movie of the resulting
% trajectories

% issues/to-do:

% short-hand for indexing coordinates
x =     1;
y =     2;
% phi =   3;

if nargin <=3
    rc = 0.035;
    display(['Setting node radius to ' num2str(rc) ' for animations.'])
end

vid = VideoWriter(filename,'MPEG-4');
open(vid)
nFrames = size(xyphiarray,4);
N = size(xyphiarray,1);
M = size(xyphiarray,2);
plotColors = lines(N);
angles = linspace(0,2*pi,8)'; % for plotting node size
nAngles = length(angles);

% set overall axes limits
xrange = minmax(reshape(xyphiarray(:,:,x,:),1,nnz(xyphiarray(:,:,x,:))));
yrange = minmax(reshape(xyphiarray(:,:,y,:),1,nnz(xyphiarray(:,:,y,:))));
xrange = [floor(xrange(1)) ceil(xrange(2))];
yrange = [floor(yrange(1)) ceil(yrange(2))];

figure
for frameCtr=1:nFrames
    plot(xyphiarray(:,:,x,frameCtr)',xyphiarray(:,:,y,frameCtr)','-',...
        'Marker','.','Color','k');
    hold on
    ax = gca;
        % plot circular domain boundaries, if given scalar domain size
    if nargin>=3&&numel(L)==1
            viscircles(ax,[0 0],L,'Color',[0.5 0.5 0.5],'LineWidth',2,'EnhanceVisibility',false);
    end
    ax.XLim = xrange;
    ax.YLim = yrange;
    ax.DataAspectRatio = [1 1 1];
    for objCtr = 1:N
        patch(xyphiarray(objCtr*ones(nAngles,1),:,x,frameCtr) + rc*cos(angles(:,ones(M,1))),...
            xyphiarray(objCtr*ones(nAngles,1),:,y,frameCtr) + rc*sin(angles(:,ones(M,1))),...
            plotColors(objCtr,:),'EdgeColor',plotColors(objCtr,:))
        % patches seems to be faster than viscircles
        %         viscircles(squeeze(xyphiarray(objCtr,:,[x y],frameCtr)),rc*ones(M,1),...
        %             'Color',plotColors(objCtr,:),'EnhanceVisibility',false,...
        %             'LineWidth',1);
    end
    writeVideo(vid,getframe)
    hold off
end
close(vid)

end
