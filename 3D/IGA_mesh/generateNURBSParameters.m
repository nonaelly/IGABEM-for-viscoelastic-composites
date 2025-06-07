function nurbsParameters = generateNURBSParameters(numSegment, shapeSegment, ...
    nodeSegment, circleCenters, p, isClosed)
% Function: generateNURBSParameters
% Description: Generates control points, knot vectors, and weights for
%   NURBS curves based on boundary shape information.
%
% Input:
%   numSegment (scalar): Number of boundary segments.
%   shapeSegment (vector): Shape of each boundary segment (0 for line, 1 for circle).
%   nodeSegment (matrix): Coordinates of the boundary nodes.
%   circleCenters (matrix): Coordinates of circle centers for circular arcs.
%   p (scalar): NURBS order.
%   isClosed (logical): Indicates whether the last boundary segment is repeated.
%
% Output:
%   nurbsParameters (struct): Struct containing parameters for NURBS curves.
%       - controlPoints: Control points of the NURBS curve.
%       - knotVector: Knot vector of the NURBS curve.
%       - weights: Weights of the control points.
%       - isClosed: Indicates whether the curve is closed or open.
%       - orderNURBS: NURBS order.
%

% Initialize knot vector
knotVector = zeros(1, numSegment * p + p + 2);

% Compute knot vector
for i = 1:numSegment-1
    knotVector(2 + p * i:p + p * i + 1) = i * ones(1, p);
end
knotVector(2 + p * numSegment:p + p * numSegment + 2) = numSegment * ones(1, p + 1);

% Initialize arrays for control points and weights
controlPoints = zeros(numSegment * p + 1, 2);
weights = zeros(numSegment * p + 1, 1);

j = 0;

% Loop over each boundary segment
for i = 1:numSegment-1
    switch shapeSegment(i)
        case 0
        % Handle line geometry
        [controlPoints(i * p - p + 1:i * p, :), weights(i * p - p + 1:i * p)] = ...
            getLineGeometry(nodeSegment, p, i, numSegment, isClosed);
        case 1
        % Handle circle geometry
        j = j + 1;
        [controlPoints(i * p - p + 1:i * p, :), weights(i * p - p + 1:i * p)] = ...
            getCircleGeometry(nodeSegment, circleCenters, p, i, j, numSegment, isClosed);
        otherwise
        error("Invalid boundary shape.");
    end
end

% Additional code for handling the last boundary segment
i = numSegment;
switch shapeSegment(i)
    case 0
    % Handle line geometry
    [controlPoints(i * p - p + 1:i * p + 1, :), weights(i * p - p + 1:i * p + 1)] = ...
        getLineGeometry(nodeSegment, p + 1, i, numSegment, isClosed);
    case 1
    % Handle circle geometry
    j = j + 1;
    [controlPoints(i * p - p + 1:i * p + 1, :), weights(i * p - p + 1:i * p + 1)] = ...
        getCircleGeometry(nodeSegment, circleCenters, p + 1, i, j, numSegment, isClosed);
    otherwise
    error("Invalid boundary shape.");
end

% Create struct to store NURBS parameters
nurbsParameters = struct;
nurbsParameters.controlPoints = controlPoints;
nurbsParameters.knotVector = knotVector;
nurbsParameters.weights = weights;
nurbsParameters.isClosed = isClosed;
nurbsParameters.orderNURBS = p;

end

%% Sub functions.
% =========================================================================
% Function to handle line geometry
function [controlPts, weights] = getLineGeometry(nodeSegment, numControl, ...
    idx, numSegment, isClosed)
controlPts = zeros(numControl, 2);
X1 = nodeSegment(idx, 1);
Y1 = nodeSegment(idx, 2);

if isClosed && numSegment == idx
    X2 = nodeSegment(1, 1);
    Y2 = nodeSegment(1, 2);
else
    X2 = nodeSegment(idx + 1, 1);
    Y2 = nodeSegment(idx + 1, 2);
end

if numSegment == idx
    p = numControl - 1;
else
    p = numControl;
end

A = (X2 - X1) / p;
B = (Y2 - Y1) / p;
for index = 1:numControl
    controlPts(index, :) = [X1 + A * index - A, Y1 + B * index - B];
end
weights = ones(numControl, 1);

end

% =========================================================================
% Function to handle circle geometry
function [controlPts, weights] = getCircleGeometry(nodeSegment, circleCenters, ...
    numControl, idx, j, numSegment, isClosed)
controlPts = zeros(numControl, 2);
X1 = nodeSegment(idx, 1);
Y1 = nodeSegment(idx, 2);

if isClosed && numSegment == idx
    X2 = nodeSegment(1, 1);
    Y2 = nodeSegment(1, 2);
else
    X2 = nodeSegment(idx + 1, 1);
    Y2 = nodeSegment(idx + 1, 2);
end
X0 = circleCenters(j, 1);
Y0 = circleCenters(j, 2);
A1 = [X1 - X0, Y1 - Y0];
A2 = [X2 - X0, Y2 - Y0];
cos2theta = dot(A1, A2) / (norm(A1) * norm(A2));
c = sqrt((1 + cos2theta) / 2);

A = (X2 - X0) * (Y1 - Y0) - (X1 - X0) * (Y2 - Y0);
X3 = -(X1 * (X1 - X0) * (Y2 - Y0) - X2 * (X2 - X0) * (Y1 - Y0) + (Y1 - Y0) * (Y2 - Y0) * (Y1 - Y2)) / A;
Y3 = -(Y2 * (X1 - X0) * (Y2 - Y0) - Y1 * (X2 - X0) * (Y1 - Y0) + (X1 - X0) * (X2 - X0) * (X2 - X1)) / A;

controlPts(1, :) = [X1, Y1];
controlPts(2, :) = [X3, Y3];

if numSegment == idx
    controlPts(3, :) = [X2, Y2];
    weights = [1; c; 1];
else
    weights = [1; c];
end

end
