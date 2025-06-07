function [jacobian, normals, r, dr, drdn] = getKernelParameters3D(element_coordinates, ...
    collocation_coordinates, N, dNu, dNv)
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
dxydxi1 = pagemtimes(dNu, element_coordinates); 
dxydxi2 = pagemtimes(dNv, element_coordinates); 

% Compute normals
normals = cross(dxydxi1,dxydxi2, 2);

jacobian = vecnorm(normals, 2, 2); 

normals = normals./jacobian;

% Compute relative distances
fieldPt = pagemtimes(N, element_coordinates); 
relDist = fieldPt - collocation_coordinates;

% Compute distance and its derivative
r = vecnorm(relDist, 2, 2);
dr = 1 ./ r .* relDist;

% Compute dot product of dr and normals
drdn = sum(dr .* normals, 2);

end
