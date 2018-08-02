function [ ] = animateTrajectories(xyphiarray,filename)
% takes in an array of N by x y phi by T and makes a movie of the resulting
% trajectories

% issues/to-do:
% - object size is not to scale

vid = VideoWriter(filename,'MPEG-4');
open(vid)
nFrames = size(xyphiarray,3);

% set overall axes limits
xrange = round(minmax(reshape(xyphiarray(:,1,:),1,nnz(xyphiarray(:,1,:)))));
yrange = round(minmax(reshape(xyphiarray(:,2,:),1,nnz(xyphiarray(:,2,:)))));

for frameCtr=1:nFrames
    quiver(xyphiarray(:,1,frameCtr),xyphiarray(:,2,frameCtr),...
        cos(xyphiarray(:,3,frameCtr)),sin(xyphiarray(:,3,frameCtr)),0.05,...% arrow scaling factor
        'Marker','o','ShowArrowHead','off');
    xlim(xrange), ylim(yrange)
    writeVideo(vid,getframe)
end
close(vid)

end

