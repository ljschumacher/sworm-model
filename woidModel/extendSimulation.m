function [] = extendSimulation(filename,Textend)
% load and resume a simulation
load(filename)
if ~exist('T','var')
    warning('Value of T unknown, using estimate.')
    T = round(size(xyarray,4)*dT*saveEvery)
    % this will throw an error if dT or saveEvery are not existing
    % variables, which could happen if they are in the param struct
end
[xyarray2, currentState] = runWoids(Textend,N,M,L,param,'resumeState',currentState);
xyarray = cat(4,xyarray,xyarray2(:,:,:,2:end));
filename = [strrep(filename,'.mat',''),'_extended.mat'];
T = T + Textend;
clear xyarray2
save(filename)
end

