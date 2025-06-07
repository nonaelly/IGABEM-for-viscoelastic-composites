function [tranU, tranT] = transformMatrix3D(nurbsStr,collocNew)
% Function: transformMatrix
% Description: Forms the transformation matrices for displacement and traction.
%
% Input:
%   nurbsStr: Structure containing boundary information.
%
% Output:
%   tranU: Transformation matrix for displacement.
%   tranT: Transformation matrix for traction.

D = 3;
tolerance = 1e-14;

% Initialize indices and counters
[indU, indT, indNumU, indNumT, indNumB] = initializeIndices(nurbsStr,collocNew);

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
    ip = indNumB(indU(i, 1)) + indU(i, 2);
    c = indU(i, 3);
    xi = nurbsStr(ip).collocXi(c, :);
    row = i * D - D + 1:i * D;

    for e = 1:size(nurbsStr(ip).elemRange,1)
        range = nurbsStr(ip).elemRange(e, :);
        uConn = nurbsStr(ip).uConn(e, :);
        % Global index on the boundary.
        uConnBou = nurbsStr(ip).globU(uConn);

        if checkPointInElement(xi, range)
            [xi_u, xi_v] = getXiParamter(xi(1), xi(2), range);
            [subDisplacement, column] = computeSubMatrix(nurbsStr(ip), [xi_u, xi_v], ...
                D, uConnBou, indNumU(indU(i, 1)));
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
    ip = indNumB(indT(i, 1)) + indT(i, 2);
    c = indT(i, 3);
    xi = nurbsStr(ip).collocXi(c, :);
    row = i * D - D + 1:i * D;

    for e = 1:size(nurbsStr(ip).elemRange,1)
        range = nurbsStr(ip).elemRange(e, :);
        tConn = nurbsStr(ip).tConn(e, :);
        % Global index on the boundary.
        tConnBou = nurbsStr(ip).globT(tConn);

        if checkPointInElement(xi, range)
            [xi_u, xi_v] = getXiParamter(xi(1), xi(2), range);
            [subTraction, column] = computeSubMatrix(nurbsStr(ip), ...
                [xi_u, xi_v], D, tConnBou, indNumT(indT(i, 1)));
            tranT(row, column) = subTraction;
            break;
            % The upper method is also acceptable.
            % 2024.12.30, I find the lower one is unacceptable. For some
            % special cases when the point is on the edge, the transform
            % matrix is not right.
%             tranT(row, column) = tranT(row, column) + subTraction;

        end
    end
end
end

%======================================================
function [subMatrix, column] = computeSubMatrix(str, xi, D, connection, offset)
% Computes the sub-matrix for displacement or traction transformation
p = str.p;
q = str.q;
weights = str.weights;
knotU = str.knotU;
knotV = str.knotV;
[N, ~, ~] = NURBS2DBasisDers(xi, p, q, knotU, knotV, weights');

% Build the sub-matrix
subMatrix = zeros(D, 3 * (p + 1) * (q + 1));
for j = 1:D
    subMatrix(j, j:D:end) = N;
end
connection = connection + offset;
column = zeros(1, 3 * (p + 1) * (q + 1));
for k = 1:D
    column(k:D:end) = connection * D - D + k;
end
end
