function [ ] = animateChainTrajectories(xyphiarray,filename)
% takes in an array of N by M by x y phi by T and makes a movie of the resulting
% trajectories

% issues/to-do:
% - node size is not to scale
% - plot segments of correct length

vid = VideoWriter(filename,'MPEG-4');
open(vid)
nFrames = size(xyphiarray,4);
% set overall axes limits
xrange = minmax(reshape(xyphiarray(:,:,1,:),1,nnz(xyphiarray(:,:,1,:))));
yrange = minmax(reshape(xyphiarray(:,:,2,:),1,nnz(xyphiarray(:,:,2,:))));
xrange = [floor(xrange(1)) ceil(xrange(2))];
yrange = [floor(yrange(1)) ceil(yrange(2))];

for frameCtr=1:nFrames
    quiver(xyphiarray(:,:,1,frameCtr),xyphiarray(:,:,2,frameCtr),...
        cos(xyphiarray(:,:,3,frameCtr)),sin(xyphiarray(:,:,3,frameCtr)),0.05,...% scaling factor, works as segment between nodes if quiver speed is 1
        'Marker','o','ShowArrowHead','off');
    xlim(xrange), ylim(yrange)
    writeVideo(vid,getframe)
end
close(vid)

end

