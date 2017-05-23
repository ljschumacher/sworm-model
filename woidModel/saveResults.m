function [ ] = saveResults( filename, savestruct )
% saves results encapsulated in a function, so that it can be run inside a
% parfor loop
% nVars = length(varargin);
% savevar = struct;
% for varCtr = 1:nVars
%     savevar.(inputname(varCtr+1)) = varargin{varCtr};
% end
save(filename,'-struct','savestruct')
end

