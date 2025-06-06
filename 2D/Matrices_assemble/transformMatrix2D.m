function [tranU, tranT] = transformMatrix2D(nurbsStr)
% Function: transformMatrix
% Description: Forms the transformation matrices for displacement and traction.
%
% Input:
%   nurbsStr: Structure containing boundary information.
%
% Output:
%   tranU: Transformation matrix for displacement.
%   tranT: Transformation matrix for traction.

D = 2;
tolerance = 1e-10;

% Initialize indices and counters
[indU, indT, indNumU, indNumT, indNumB] = initializeIndices2D(nurbsStr);

% Construct the displacement transformation matrix
tranU = transfromMatU(nurbsStr, indU, indNumU, indNumB, D);

% Construct the traction transformation matrix
tranT = transfromMatT(nurbsStr, indT, indNumT, indNumB, D);

% Set small values to zero
tranU(abs(tranU) < tolerance) = 0;
tranT(abs(tranT) < tolerance) = 0;
end

%======================================================
function tranU = transfromMatU(nurbsStr, indU, indNumU, indNumB, D)
% Initialize transformation matrices
numU = indNumU(end);
tranU = zeros(numU * D);
for i = 1:numU
    ip = indU(i, 1);
    c = indU(i, 2);
    xi = nurbsStr(ip).collocXi(c, :);
    row = i * D - D + 1:i * D;

    for e = 1:size(nurbsStr(ip).elemRange,1)
        range = nurbsStr(ip).elemRange(e, :);
        uConn = nurbsStr(ip).uConn(e, :);
        gloIdx = nurbsStr(ip).globIdx(e, :);

        % Judge the type of integral
        [xi_u, flag] = JudgeSingular2D(nurbsStr, xi, range, ip, ip, e);

        if flag
            [subDisplacement, column] = computeSubMatrix(nurbsStr(ip), gloIdx, ...
                xi_u, D, uConn, indNumU(indU(i, 1)));
            tranU(row, column) = subDisplacement;
            break;
        end
    end
end
end

%======================================================
function tranT = transfromMatT(nurbsStr, indT, indNumT, indNumB, D)
numT = indNumT(end);
tranT = zeros(numT * D);
for i = 1:numT
    ip = indT(i, 1);
    t2u = nurbsStr(ip).t2uConn;
    c = t2u(indT(i, 2));
    xi = nurbsStr(ip).collocXi(c, :);
    row = i * D - D + 1:i * D;

    for e = 1:size(nurbsStr(ip).elemRange,1)
        range = nurbsStr(ip).elemRange(e, :);
        tConn = nurbsStr(ip).tConn(e, :);
        gloIdx = nurbsStr(ip).globIdx(e, :);

        % Judge the type of integral
        [xi_u, flag] = JudgeSingular2D(nurbsStr, xi, range, ip, ip, e);

        flagT = (ismember(i,tConn+indNumT(indT(i, 1)))) && flag;

        if flagT
            [subTraction, column] = computeSubMatrix(nurbsStr(ip), gloIdx, ...
                xi_u, D, tConn, indNumT(indT(i, 1)));
            tranT(row, column) = subTraction;
            break;
        end
    end
end
end

%======================================================
function [subMatrix, column] = computeSubMatrix(str, gloIdx, xi, D, connection, offset)
% Computes the sub-matrix for displacement or traction transformation
p = str.p;
weights = str.weights;
knotU = str.knotU;

[N, ~] = NURBSbasis(gloIdx, p, xi, knotU, weights');

% Build the sub-matrix
subMatrix = zeros(D, D * (p + 1));
for j = 1:D
    subMatrix(j, j:D:end) = N;
end
connection = connection + offset;
column = zeros(1, D * (p + 1));
for k = 1:D
    column(k:D:end) = connection * D - D + k;
end
end
