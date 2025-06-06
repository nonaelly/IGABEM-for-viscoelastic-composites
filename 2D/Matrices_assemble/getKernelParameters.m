function [jacobian, normals, r, dr, drdn] = getKernelParameters(element_coordinates, collocation_coordinates, N, dN)
% getKernelParameters: Calculates the parameters needed to evaluate the kernels.
%
% Inputs:
%   - element_coordinates: Coordinates of the element.
%   - collocation_coordinates: Coordinates of the collocation points.
%   - N: NURBS basis functions.
%   - dN: Derivatives of NURBS basis functions.
%
% Outputs:
%   - jacobian: Jacobian of the transformation.
%   - normals: Normal vectors.
%   - r: Distance between field points and collocation coordinates.
%   - dr: Derivative of r.
%   - drdn: Dot product of dr and normals.

% =========================================================================
% Compute geometry derivatives
dxydxi = pagemtimes(dN, element_coordinates); 
jacobian = vecnorm(dxydxi, 2, 2); % Norm in every row

% Compute normals
normals = [dxydxi(:, 2, :), - dxydxi(:, 1, :)] ./ jacobian;

% Compute relative distances
fieldPt = pagemtimes(N, element_coordinates); 
relDist = fieldPt - collocation_coordinates;

% Compute distance and its derivative
r = vecnorm(relDist, 2, 2);
dr = 1 ./ r .* relDist;

% Compute dot product of dr and normals
drdn = sum(dr .* normals, 2);

end
