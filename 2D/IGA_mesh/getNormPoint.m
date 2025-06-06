function [normNode, jacobian] = getNormPoint(coordElemNodes, xi, basisIdx, ...
    knotVec, weights, p)
% Function: getNormPoint
% Description: Calculates the normal vector at a local coordinate on an element.
%
% Inputs:
%   coordElemNodes (matrix): Coordinates of the element nodes.
%   xi (scalar): Local coordinate on the element.
%   basisIdx (vector): Global basis index of the NURBS basis functions.
%   knotVec (vector): Knot vector of the NURBS curve.
%   weights (vector): Weights of the control points.
%   p (scalar): Degree of the NURBS curve.
%
% Outputs:
%   normNode (vector): Normal vector at the given local coordinate.
%   jacobian (scalar): Jacobian determinant at the given local coordinate.

% Calculate NURBS basis functions and their derivatives at xi
[~, dN] = NURBSbasis(basisIdx, p, xi, knotVec, weights);

% Calculate the geometry derivatives
dxydxi = dN * coordElemNodes;

% Calculate the Jacobian determinant
jacobian = norm(dxydxi);

% Calculate the normal vector
normNode = 1 / jacobian * [dxydxi(2), -dxydxi(1)];

end
