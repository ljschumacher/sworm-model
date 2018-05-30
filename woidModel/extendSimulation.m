function [] = extendSimulation(thisfilename,Textend,makeMov)
if nargin<3
    makeMov=false;
end
% load and resume a simulation
load(thisfilename)
if ~exist('T','var')
    warning('Value of T unknown, using estimate.')
    T = round(size(xyarray,4)*param.dT*param.saveEvery)
end
if ~exist('bc','var')
    warning('Boundary conditions unspecified, using periodic')
    bc = 'periodic';
end
[xyarray2, currentState, food2] = runWoids(Textend,N,M,L,param,...
    'bc',bc,'resumeState',currentState);
xyarray = cat(4,xyarray,xyarray2(:,:,:,2:end));
if exist('food','var')
    food = cat(3,food,food2(:,:,2:end));
end
thisfilename = [strrep(thisfilename,'.mat',''),'_extended.mat'];
T = T + Textend;
clear xyarray2
clear food2
clear Textend
save(thisfilename)

if makeMov
    addpath('visualisation')
    animateWoidTrajectories(xyarray,strrep(thisfilename,'.mat',''),L,0.035,food);
end
end

