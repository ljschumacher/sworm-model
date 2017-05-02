function [ ] = animateWoidTrajectories(xyarray,filename,L,rc)
% takes in an array of N by M by x y  by T and makes a movie of the resulting
% trajectories

% issues/to-do:

% short-hand for indexing coordinates
x =     1;
y =     2;

if nargin <=3
    rc = 0.035;
    display(['Setting node radius to ' num2str(rc) ' for animations.'])
end

if ~ismac
    vid = VideoWriter(filename,'Motion JPEG AVI');
else
    vid = VideoWriter(filename,'MPEG-4');
end
open(vid)
nFrames = size(xyarray,4);
N = size(xyarray,1);
M = size(xyarray,2);
plotColors = lines(N);
angles = linspace(0,2*pi,10)'; % for plotting node size

% set overall axes limits
xrange = minmax(reshape(xyarray(:,:,x,:),1,numel(xyarray(:,:,x,:))));
yrange = minmax(reshape(xyarray(:,:,y,:),1,numel(xyarray(:,:,y,:))));
xrange = [floor(xrange(1)) ceil(xrange(2))];
yrange = [floor(yrange(1)) ceil(yrange(2))];

figure
for frameCtr=1:nFrames
    if M>1
        plot(xyarray(:,:,x,frameCtr)',xyarray(:,:,y,frameCtr)','-',...
            'Marker','.','Color','k');
    else
        plot(xyarray(:,:,x,frameCtr)',xyarray(:,:,y,frameCtr)','.','Color','k');
    end
    hold on
    if N==1 % plot tracks for single worm
        plot(squeeze(xyarray(:,:,x,1:frameCtr)),squeeze(xyarray(:,:,y,1:frameCtr)),'-',...
            'Color',[0.5 0.5 0.5]);
    end
    ax = gca;
    % plot circular domain boundaries, if given scalar domain size
    if nargin>=3&&numel(L)==1
        viscircles(ax,[0 0],L,'Color',[0.5 0.5 0.5],'LineWidth',2,'EnhanceVisibility',false);
        ax.XLim = [-L L];
        ax.YLim = [-L L];
    else
        ax.XLim = xrange;
        ax.YLim = yrange;
    end
    ax.DataAspectRatio = [1 1 1];
    ax.XTick = [];
    ax.YTick = [];
    for objCtr = 1:N
        patch(xyarray(objCtr,:,x,frameCtr) + rc*cos(angles),...
            xyarray(objCtr,:,y,frameCtr) + rc*sin(angles),...
            plotColors(objCtr,:),'EdgeColor',plotColors(objCtr,:))
        % patches seems to be faster than viscircles
        %         viscircles(squeeze(xyarray(objCtr,:,[x y],frameCtr)),rc*ones(M,1),...
        %             'Color',plotColors(objCtr,:),'EnhanceVisibility',false,...
        %             'LineWidth',1);
    end
    writeVideo(vid,getframe)
    hold off
end
close(vid)

end
