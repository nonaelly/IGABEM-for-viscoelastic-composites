function [ctrlPts, knotU, knotV, weights] = hRefinement3D(p, q, refinement, knotU, knotV, ctrlPts, weights)
% hRefinement3D Refines a 3D NURBS surface by adding knots and redistributing control points.
%
% Inputs:
%   p (scalar): Degree of the NURBS surface in the U direction.
%   q (scalar): Degree of the NURBS surface in the V direction.
%   refinement (scalar): Refinement level.
%   knotU (vector): Knot vector in the U direction.
%   knotV (vector): Knot vector in the V direction.
%   ctrlPts (matrix): Control points of the NURBS surface.
%   weights (vector): Weights of the control points.
%
% Outputs:
%   newCtrlPts (matrix): Refined control points of the NURBS surface.
%   newKnotU (vector): Refined knot vector in the U direction.
%   newKnotV (vector): Refined knot vector in the V direction.
%   newWeights (vector): Refined weights of the NURBS surface.

% Normalize knot vectors
knotU = knotU / max(knotU);
knotV = knotV / max(knotV);

% Initial
m = length(knotU) - p - 1;
n = length(knotV) - q - 1;

% h-Refinement
% U-direction
mNew = (length(unique(knotU))-1)*(refinement(1) - 1)+m;
ctrlPtsNew = zeros(mNew*n, 3);
weightsNew = zeros(mNew*n, 1);

for j=1:n
    ind1 = m*(j-1)+1 : j*m;
    ind2 = mNew*(j-1)+1 : j*mNew;
    [ctrlPtsNew(ind2,:), knotUNew, weightsNew(ind2)] = hRefinement(p, refinement(1), ...
    knotU, ctrlPts(ind1,:), weights(ind1,:));
end
weights = weightsNew;
ctrlPts = ctrlPtsNew;
knotU = knotUNew;

% V-direction
nNew = (length(unique(knotV))-1)*(refinement(2) - 1)+n;
ctrlPtsNew = zeros(mNew*nNew, 3);
weightsNew = zeros(mNew*nNew, 1);

for j=1:mNew
    ind1 = j : mNew : n*mNew;
    ind2 = j : mNew : nNew*mNew;
    [ctrlPtsNew(ind2,:), knotVNew, weightsNew(ind2)] = hRefinement(q, refinement(2), ...
    knotV, ctrlPts(ind1,:), weights(ind1,:));
end
weights = weightsNew;
ctrlPts = ctrlPtsNew;
knotV = knotVNew;

end


