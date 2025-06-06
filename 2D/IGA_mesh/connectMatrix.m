function [globIdx, dispConn, tracConn, elemRange] = connectMatrix(p, knotVec, ctrlPts, weights, isClosed)
% Function: connectMatrix
% Description: Get the connect matrix for the global basis index (a),
% displacement (U) and traction (T)
%
% Input:
%   p (scalar): Degree of the NURBS curve.
%   knotVec (vector): Knot vector of the NURBS curve.
%   ctrlPts (matrix): Control points of the NURBS curve.
%   weights (vector): Weights of the control points.
%   isClosed (logical): Indicates whether the curve is closed or open.
%
% Output:
%   globIdx (matrix): Connect matrix for global basis functions.
%   dispConn (matrix): Connect matrix for displacement.
%   tracConn (matrix): Connect matrix for displacement.
%   elemRange (matrix): Matrix for element range.
%
% -------------------------------------------------------------------
%   element index (e)                   global basis index (a) for b1, b2, b3
%   1                                   1, 2, 3
%   2                                   ...
%   ...                                 ...
%   number of element (numElem)         ...

% Initialize
uniqueKnots = unique(knotVec);
numElem = length(uniqueKnots) - 1; % Number of elements

elemRange = zeros(numElem, 2);        % Initialize matrices
dispConn = zeros(numElem, p + 1);
elemKnotIdx = zeros(numElem, 2);
tracConn = zeros(numElem, p + 1);

% Determine element ranges and the corresponding knot indices
e = 1;
previousKnotVal = 0;
tol = 1e-8;

for i = 1:length(knotVec)
    currentKnotVal = knotVec(i);
    if knotVec(i) ~= previousKnotVal
        elemRange(e, :) = [previousKnotVal, currentKnotVal];
        elemKnotIdx(e, :) = [i - 1, i];
        e = e + 1;
    end
    previousKnotVal = currentKnotVal;
end

numRepeatedKnots = 0;

for e = 1:numElem
    indices = (elemKnotIdx(e, 1) - p + 1) : elemKnotIdx(e, 1);
    previousKnotVals = knotVec(indices);
    currentKnotVals = ones(1, p) * knotVec(elemKnotIdx(e, 1));
    if isequal(previousKnotVals, currentKnotVals) && length(nonzeros(previousKnotVals)) > 1
        % Judge the norm vector of both sides.
        % Left
        xi = knotVec(elemKnotIdx(e, 1));
        idxLeft = (elemKnotIdx(e - 1, 1) - p) : elemKnotIdx(e - 1, 1);
        [~, dNLeft] = NURBSbasis(idxLeft, p, xi - eps, knotVec, weights);
        % Tangent vector 
        tanLeft = dNLeft * ctrlPts(idxLeft, :);
        normLeft = [tanLeft(2), -tanLeft(1)] / norm(tanLeft);

        % Right
        idxRight = (elemKnotIdx(e, 1) - p) : elemKnotIdx(e, 1);
        [~, dNRight] = NURBSbasis(idxRight, p, xi + eps, knotVec, weights);
        tanRight = dNRight * ctrlPts(idxRight, :);
        normRight = [tanRight(2), -tanRight(1)] / norm(tanRight);
        if norm(normLeft - normRight) > tol
            numRepeatedKnots = numRepeatedKnots + 1;
        end
    end
    dispConn(e, :) = (elemKnotIdx(e, 1) - p) : elemKnotIdx(e, 1);
    tracConn(e, :) = dispConn(e, :) + numRepeatedKnots;
end

globIdx = dispConn;

if isClosed
    dispConn(end, end) = 1;  % the last point is equal to the first point
    % Smooth at last point
    % Left
    idxLeft = globIdx(end, :);
    [~, dNLeft] = NURBSbasis(idxLeft, p, 1 - eps, knotVec, weights);
    % Tangent vector
    tanLeft = dNLeft * ctrlPts(idxLeft, :);
    normLeft = [tanLeft(2), -tanLeft(1)] / norm(tanLeft);

    % Right
    idxRight = globIdx(1, :);
    [~, dNRight] = NURBSbasis(idxRight, p, 0 + eps, knotVec, weights);
    tanRight = dNRight * ctrlPts(idxRight, :);
    normRight = [tanRight(2), -tanRight(1)] / norm(tanRight);
    if norm(normLeft - normRight) < tol
        tracConn(end, end) = 1;
    end
end

end
