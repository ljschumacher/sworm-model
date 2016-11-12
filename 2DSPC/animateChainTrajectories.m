function [ ] = animateChainTrajectories(xyphiarray,filename)
% takes in an array of N by M by x y phi by T and makes a movie of the resulting
% trajectories

% issues/to-do:
% - node size is not to scale

vid = VideoWriter(filename,'MPEG-4');
open(vid)
nFrames = size(xyphiarray,4);
% set overall axes limits
xrange = minmax(reshape(xyphiarray(:,:,1,:),1,nnz(xyphiarray(:,:,1,:))));
yrange = minmax(reshape(xyphiarray(:,:,2,:),1,nnz(xyphiarray(:,:,2,:))));
xrange = [floor(xrange(1)) ceil(xrange(2))];
yrange = [floor(yrange(1)) ceil(yrange(2))];

for frameCtr=1:nFrames
    plot(xyphiarray(:,:,1,frameCtr)',xyphiarray(:,:,2,frameCtr)','-',...
        'Marker','.','Color','k');
    ax = gca;
    ax.XLim = xrange;
    ax.YLim = yrange;
    ax.DataAspectRatio = [1 1 1];
    writeVideo(vid,getframe)
end
close(vid)

end

