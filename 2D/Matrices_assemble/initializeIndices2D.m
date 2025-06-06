function [indU, indT, indNumU, indNumT, indNumB] = initializeIndices2D(nurbsStr)
% Initialize index matrices for displacement and traction
%
% Inputs:
%   nurbsStr: Structure array containing boundary element data
%
% Outputs:
%   indU: Index vector for displacement
%   indT: Index vector for traction
%   indNumU: Index vector for number of displacement on different boundary.
%   indNumT: Index vector for number of traction on different boundary.
%   indNumB: Index vector for number of boundary

numBoundaries = length(nurbsStr);

% Initialize cumulative counters
indNumU = zeros(numBoundaries + 1, 1);
indNumT = zeros(numBoundaries + 1, 1);
indNumB = zeros(numBoundaries + 1, 1);

% Preallocate index matrices
indU = [];
indT = [];

% Iterate through each boundary structure
for i = 1:numBoundaries
    numCollocU = nurbsStr(i).numCollocU;
    numCollocT = nurbsStr(i).numCollocT;

    % Update cumulative counters
    indNumU(i+1) = indNumU(i) + numCollocU;
    indNumT(i+1) = indNumT(i) + numCollocT;
    indNumB(i+1) = indNumB(i) + 1;
    
    % Create index vector for displacement and traction
    indU = [indU; [i*ones(numCollocU,1), (1:numCollocU)']];
    indT = [indT; [i*ones(numCollocT,1), (1:numCollocT)']];
end

end
