function [ ] = animateWoidTrajectories(xyarray,filename,L,rc)
% takes in an array of N by M by x y  by T and makes a movie of the resulting
% trajectories

% issues/to-do:

% short-hand for indexing coordinates
x =     1;
y =     2;

if nargin <4
    rc = 0.035;
    disp(['Setting node radius to ' num2str(rc) ' for animations.'])
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
angles = linspace(0,2*pi,20)'; % for plotting node size

% set overall axes limits
xrange = minmax(reshape(xyarray(:,:,x,:),1,numel(xyarray(:,:,x,:))));
yrange = minmax(reshape(xyarray(:,:,y,:),1,numel(xyarray(:,:,y,:))));
% xrange = [floor(xrange(1)) ceil(xrange(2))];
% yrange = [floor(yrange(1)) ceil(yrange(2))];

% calculate markerSize to plot (in points) from
figure
for frameCtr=1:nFrames
    if M>1% plot connecting lines btw nodes
        % don't plot connecting lines for objects that span across a
        % periodic boundary
        excludedObjects = any(any(abs(diff(xyarray(:,:,:,frameCtr),1,2))>min(L./2),3),2);
        plot(xyarray(~excludedObjects,:,x,frameCtr)',xyarray(~excludedObjects,:,y,frameCtr)','k-');
        if any(~excludedObjects), hold on, end
    elseif rc==0
        plot(xyarray(:,:,x,frameCtr)',xyarray(:,:,y,frameCtr)','.','Color','k');
        hold on
    end
    if N==1 % plot tracks for single worm
        plot(squeeze(xyarray(:,:,x,1:frameCtr)),squeeze(xyarray(:,:,y,1:frameCtr)),'-',...
            'Color',[0.5 0.5 0.5]);
        hold on
    end
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
    if rc>0
        for objCtr = 1:N
            patch(xyarray(objCtr,:,x,frameCtr) + rc*cos(angles),...
                xyarray(objCtr,:,y,frameCtr) + rc*sin(angles),...
                plotColors(objCtr,:),'EdgeColor','none')
            if objCtr==1
                hold on
            end
            % patches seems to be faster than viscircles
            %         viscircles(squeeze(xyarray(objCtr,:,[x y],frameCtr)),rc*ones(M,1),...
            %             'Color',plotColors(objCtr,:),'EnhanceVisibility',false,...
            %             'LineWidth',1);
        end
    end
    % plot heads
    plot(xyarray(:,1,x,frameCtr)',xyarray(:,1,y,frameCtr)','.','Color','k');
    writeVideo(vid,getframe)
    hold off
end
close(vid)

end