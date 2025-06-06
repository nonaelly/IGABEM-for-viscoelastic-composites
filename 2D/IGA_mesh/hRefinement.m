function [controlPoints, knotVector, weights] = hRefinement(p, refinement, ...
    knotVector, controlPoints, weights)
% Function: hRefinement
% Description: Refines the 1D NURBS curve by adding knots and redistributing 
%   control points.
% Reference: The NURBS Book (Page 104-105).
%
% Input:
%   p (scalar): Degree of the NURBS curve.
%   refinement (scalar): Number of refinements to perform.
%   knotVector (vector): Knot vector of the NURBS curve.
%   controlPoints (matrix): Control points of the NURBS curve.
%   weights (vector): Weights of the control points.
%
% Output:
%   controlPoints (matrix): Refined control points of the NURBS curve.
%   knotVector (vector): Refined knot vector of the NURBS curve.
%   weights (vector): Refined weights of the NURBS curve.
%

% Dimension
d = size(controlPoints,2);

% Normalize knot vector
normalizedKnotVector = knotVector / max(knotVector);

% Weighted control points. (wi*xi,wi*yi,wi*zi,wi)
weightedControlPoints = [controlPoints .* weights, weights];

% Identify unique knots
uniqueKnots = unique(normalizedKnotVector);

% Initialize matrix to store new knots
if length(refinement) ~= length(uniqueKnots) - 1
    error('Refinement not match')
end
newKnots = zeros(sum(refinement - 1), 1);

% Generate new knots
temp = 0;
for i = 1 : length(uniqueKnots) - 1
    distribution = linspace(uniqueKnots(i), uniqueKnots(i + 1), refinement(i) + 1);
    ind = 1 + temp : refinement(i) - 1 + temp;
    newKnots(ind, :) = distribution(2:end-1);
    temp = temp + refinement(i) - 1;
end

% Reshape and sort new knots
newKnots = sort(reshape(newKnots, 1, []));

oldKnotVector = normalizedKnotVector;

% Perform refinement
for j = 1:length(newKnots)
    % Initialize new control points array
    newControlPoints = zeros(size(weightedControlPoints, 1) + 1, d+1);
    uBar = newKnots(j);
    k = find(oldKnotVector > uBar, 1) - 1;
    for i = 1:size(newControlPoints, 1)
        if i <= (k - p)
            alpha = 1;
        elseif i >= (k - p + 1) && i <= k
            alpha = (uBar - oldKnotVector(i)) / (oldKnotVector(i + p) - oldKnotVector(i));
        else
            alpha = 0;
        end

        newPoint = zeros(1, d+1);
        if i ~= 1
            newPoint = (1 - alpha) * weightedControlPoints(i - 1, :);
        end
        if i ~= length(newControlPoints)
            newPoint = newPoint + alpha * weightedControlPoints(i, :);
        end

        newControlPoints(i, :) = newPoint;
    end
    weightedControlPoints = newControlPoints;
    knotVector = sort([knotVector, uBar]);
    oldKnotVector = knotVector;
end

% Divide weights from control points
weights = weightedControlPoints(:, d+1);
controlPoints = weightedControlPoints(:, 1:d) ./ weights;

end
