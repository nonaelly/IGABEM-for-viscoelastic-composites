function [globIdx, uConn, tConn, t2uConn, elemRange] = connectMatrix3D(p, q, knotU, knotV, isClosedU, isClosedV)
% connectMatrix3D Generates connectivity matrices for 3D NURBS surface.
%
% Inputs:
%   p (scalar): Degree of the NURBS surface in the U direction.
%   q (scalar): Degree of the NURBS surface in the V direction.
%   knotU (vector): Knot vector in the U direction.
%   knotV (vector): Knot vector in the V direction.
%   isClosedU (logical): Indicates whether the U direction is closed.
%   isClosedV (logical): Indicates whether the V direction is closed.
%
% Outputs:
%   globIdx (matrix): Global basis index matrix.
%   uConn (matrix): Displacement connection matrix.
%   tConn (matrix): Traction connection matrix.
%   elemRange (matrix): Element range matrix.

% -------------------------------------------------------------------------
% Generate U direction connection matrix
[globIdxU, uConnU, tConnU, elemRangeU] = generateConnectionMatrix(p, knotU, isClosedU);

% Generate V direction connection matrix
[globIdxV, uConnV, tConnV, elemRangeV] = generateConnectionMatrix(q, knotV, isClosedV);

% Combine U and V direction connection matrices
numElemU = size(globIdxU, 1);
numElemV = size(globIdxV, 1);
m = length(knotU) - p - 1;
n = length(knotV) - q - 1;

matrixGlobal = reshape(1:m*n, m, n)';
uGlobal = (1:max(uConnU(end,:))) + max(uConnU(end,:))*((1:max(uConnV(end,:))) - 1)';
tGlobal = (1:max(tConnU(end,:))) + max(tConnU(end,:))*((1:max(tConnV(end,:))) - 1)';

globIdx = zeros(numElemU*numElemV, (p+1)*(q+1));
uConn = zeros(numElemU*numElemV, (p+1)*(q+1));
tConn = zeros(numElemU*numElemV, (p+1)*(q+1));
elemRange = zeros(numElemU*numElemV, 4);
t2uConn = zeros(max(tConn(:)), 1);

for i=1:numElemV
    for j=1:numElemU
        ind = j + (i-1)*numElemU;
        globIdx(ind, :) = reshape((matrixGlobal(globIdxV(i,:), globIdxU(j,:)))', 1, []);
        uConn(ind, :) = reshape(uGlobal(uConnV(i,:), uConnU(j,:))', 1, []);
        tConn(ind, :) = reshape(tGlobal(tConnV(i,:), tConnU(j,:))', 1, []);
        elemRange(ind, :) = [elemRangeU(j, :), elemRangeV(i, :)];
        t2uConn(tConn(ind,:)) = uConn(ind,:);
    end
end

end

%======================================================
function [globIdx, uConn, tConn, elemRange] = generateConnectionMatrix(p, knot_vector, isClosed)
% generateConnectionMatrix Generates the connection matrix for a given knot vector.
%
% Inputs:
%   p (scalar): Degree of the NURBS curve.
%   knot_vector (vector): Knot vector of the NURBS curve.
%   isClosed (logical): Indicates whether the curve is closed or open.
%
% Outputs:
%   globIdx (matrix): Global basis index matrix.
%   uConn (matrix): Displacement connection matrix.
%   tConn (matrix): Traction connection matrix.
%   elemRange (matrix): Element range matrix.

% Initialize variables
uniqueKnots = unique(knot_vector);
num_elements = length(uniqueKnots) - 1;

elemRange = zeros(num_elements, 2);  % Initialize element range matrix
uConn = zeros(num_elements, p+1);
knotIndices = zeros(num_elements, 2);
tConn = zeros(num_elements, p+1);

% Determine element ranges and the corresponding knot indices
e = 1;
previousKnotVal = 0;

for i = 1:length(knot_vector)
    currentKnotVal = knot_vector(i);
    if knot_vector(i) ~= previousKnotVal
        elemRange(e, :) = [previousKnotVal, currentKnotVal];
        knotIndices(e, :) = [i-1, i];
        e = e + 1;
    end
    previousKnotVal = currentKnotVal;
end

numRepeatedKnots = 0;

for e = 1:num_elements
    indices = (knotIndices(e, 1) - p + 1):knotIndices(e, 1);
    previousKnotVals = knot_vector(indices);
    currentKnotVals = ones(1, p) * knot_vector(knotIndices(e, 1));
    if isequal(previousKnotVals, currentKnotVals) && length(nonzeros(previousKnotVals)) > 1
        % Assume it is smooth
        % numRepeatedKnots = numRepeatedKnots + 1;
    end
    uConn(e, :) = (knotIndices(e, 1) - p):knotIndices(e, 1);
    tConn(e, :) = uConn(e, :) + numRepeatedKnots;
end

globIdx = uConn;

if isClosed
    uConn(end, end) = 1;  % the last point is equal to the first point
end

end
