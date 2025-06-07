function [preU, preT, unU, unT, vecT, vecU, face325U, face146U, face325T, face146T] ...
    = setPeriodicBoundary3D(nurbsStr, collocNew, epsilonPB)
% Sets the periodic boundary conditions for a 3D problem.
%
% Inputs:
%   nurbsStr: Structure array containing boundary element data.
%   collocNew: Struct containing new collocation points and global indices.
%   epsilonPB: Matrix for periodic boundary conditions.
%
% Outputs:
%   preU: Prescribed displacement degrees of freedom.
%   preT: Prescribed traction degrees of freedom.
%   unU: Unknown displacement degrees of freedom.
%   unT: Unknown traction degrees of freedom.
%   collocT: Collected traction values.
%   collocU: Collected displacement values.
%   periodU325: Periodic indices for displacement on faces 3, 2, 5.
%   periodU146: Periodic indices for displacement on faces 1, 4, 6.
%   periodT325: Periodic indices for traction on faces 3, 2, 5.
%   periodT146: Periodic indices for traction on faces 1, 4, 6.

% Initialize
D = 3;
[indU, indT, indNumU, indNumT, indNumB] = initializeIndices(nurbsStr, collocNew);
numU = indNumU(end);
numT = indNumT(end);
numBou = length(collocNew);
[face325U, face146U, face325T, face146T] = deal(cell(D,1));
preU = [];
preT = [];

rveU = indU(indU(:,1) == 1,:);
rveT = indT(indT(:,1) == 1,:);

known = [3,2,5];
unknown = [1,4,6];
for i=1:D

    % indNumB(rveU(:,1)) = 0.
    u325 = nurbsStr(known(i)).globU( rveU(rveU(:,2) == known(i), 3) );
    u146 = nurbsStr(unknown(i)).globU( rveU(rveU(:,2) == unknown(i), 3) );
    t325 = nurbsStr(known(i)).globT( rveT(rveT(:,2) == known(i), 3) );
    t146 = nurbsStr(unknown(i)).globT( rveT(rveT(:,2) == unknown(i), 3) );

    % Keep the order of opposite face same.
    switch i
        case {1,2}
            numeta = nurbsStr(i).refinement(2) + 2;
            numxi = nurbsStr(i).refinement(1) + 2;
            tempu = zeros(length(u325),1);
            tempt = zeros(length(t325),1);
            for j=1:numeta
                tempu(1+(j-1)*numxi : j*numxi) = flip(u325(1+(j-1)*numxi : j*numxi));
                tempt(1+(j-1)*numxi : j*numxi) = flip(t325(1+(j-1)*numxi : j*numxi));
            end
            u325 = tempu;
            t325 = tempt;
        case 3
            numeta = nurbsStr(i).refinement(2) + 2;
            numxi = nurbsStr(i).refinement(1) + 2;
            tempu = zeros(length(u325),1);
            tempt = zeros(length(t325),1);
            for j=1:numeta
                tempu(1+(j-1)*numxi : j*numxi) = u325(1+(numeta-j)*numxi : (numeta-j+1)*numxi);
                tempt(1+(j-1)*numxi : j*numxi) = t325(1+(numeta-j)*numxi : (numeta-j+1)*numxi);
            end
            u325 = tempu;
            t325 = tempt;
    end

    n = size(u325,1);
    [indU325,indU146,indT325,indT146] = deal(zeros(D*n,1));
    for k=1:D
        indU325(k:D:end) = u325*D-D+k;
        indU146(k:D:end) = u146*D-D+k;
        indT325(k:D:end) = t325*D-D+k;
        indT146(k:D:end) = t146*D-D+k;
    end

    face325U{i} = indU325;
    face146U{i} = indU146;
    face325T{i} = indT325;
    face146T{i} = indT146;
    preU = [preU; indU325];
    preT = [preT; indT325];
end

% Inclusion's traction.
numRVET = sum(indT(:,1) == 1);
preT = [preT; (1+numRVET*D : numT*D)'];

preU = unique(preU);
preT = unique(preT);
unU = setdiff((1:numU*D)',preU);
unT = setdiff((1:numT*D)',preT);

% Discontinuous element so NO conner points.
vecU = zeros(D*numU,1);
vecT = zeros(D*numT,1);
end