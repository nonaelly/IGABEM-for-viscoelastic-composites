function [H, G] = integralRegular3D(elCoords, collocCoords, range, p, q, ...
    knotU, knotV, weights, const1, const2, const3, const4, gptUR, gptVR, gwtUR, gwtVR, D, id)
% integralRegular3D: Computes the H and G matrices for the regular case.
%
% Inputs:
%   elCoords: Coordinates of the element.
%   collocCoords: Coordinates of the source points.
%   range: Range of the element in parameter space [xi_i, xi_i+1].
%   p, q: Degrees of the NURBS basis functions.
%   knotU, knotV: Knot vectors of the NURBS basis functions.
%   weights: Weights of the control points.
%   const1, const2, const3, const4: Material constants.
%   gptUR, gptVR: Gauss points in the parameter space.
%   gwtUR, gwtVR: Gauss weights in the parameter space.
%   D: Dimension of the problem.
%
% Outputs:
%   H, G: The H and G matrices for the regular case.

% Jacobian from parent to parameter space
jacobParam = (range(2) - range(1)) * (range(4) - range(3)) / 4;

% Gauss points in parameter space
xiU = convertToParamSpace(gptUR(:), range(1:2));
xiV = convertToParamSpace(gptVR(:), range(3:4));
numGauss = length(xiU); % Number of Gauss points

N = zeros(1, (p+1)*(q+1), numGauss); % NURBS basis functions at Gauss points
dNu = N; % Derivatives w.r.t. U
dNv = N; % Derivatives w.r.t. V

% Compute NURBS basis functions and derivatives at Gauss points
for i = 1:numGauss
    [N(:,:,i), dNu(:,:,i), dNv(:,:,i)] = NURBS2DBasisDers([xiU(i,:), xiV(i,:)], ...
        p, q, knotU, knotV, weights');
end

% Get kernel parameters
[jacobXi, normals, r, dr, drdn] = getKernelParameters3D(elCoords, collocCoords, N, dNu, dNv);

% Adjust signs for inclusion elements
if id > 6
    normals = -normals;
    drdn = -drdn;
end

% Compute final Jacobian
jacob = jacobXi * jacobParam; % The final Jacobian we use

% Compute H and G matrices using ElastKernel
[H, G] = ElastKernel(r, dr, drdn, const1, const2, const3, const4, ...
    numGauss, jacob, normals, N, gwtUR, gwtVR, D, p, q);

end

% =========================================================================
function [H, G] = ElastKernel(r, dr, drdn, const1, const2, const3, const4, ...
    numGauss, jacob, normals, N, gwtU, gwtV, D, p, q)
% ElastKernel: Calculates the elastic kernel matrix.
%
% Inputs:
%   r: Distance between field points and source coordinates.
%   dr: Derivative of r.
%   drdn: Dot product of dr and normals.
%   const1, const2, const3, const4: Constants for elasticity calculation.
%   numGauss: Number of Gaussian points.
%   jacob: Jacobian of the transformation.
%   normals: Normal vectors.
%   N: NURBS basis functions.
%   gwtU, gwtV: Gaussian weights.
%   D: Dimensionality of the problem.
%   p, q: Degrees of the NURBS basis functions.
%
% Outputs:
%   H, G: Elastic kernel matrices.

numCtrlPts = (p+1)*(q+1);

% Reshape and broadcast matrices for element-wise operations
dr = repmat(dr, [(p+1)*(q+1), 1, 1]);
dr = permute(dr, [2, 1, 3]);
drExp = reshape(dr, [D, 1, numCtrlPts, numGauss]);
drTExp = permute(drExp, [2 1 3 4]);

normals = repmat(normals, [(p+1)*(q+1), 1, 1]);
normals = permute(normals, [2, 1, 3]);
normExp = reshape(normals, [D, 1, numCtrlPts, numGauss]);

r = repmat(r, [(p+1)*(q+1), 1, 1]);
rExp = reshape(r, [1, 1, numCtrlPts, numGauss]);

drdn = repmat(drdn, [(p+1)*(q+1), 1, 1]);
drdnExp = reshape(drdn, [1, 1, numCtrlPts, numGauss]);

jacob = reshape(jacob, [1, 1, 1, numGauss]);
jacobExp = repmat(jacob, [1, 1, numCtrlPts, 1]);

NExp = reshape(N, [1, 1, numCtrlPts, numGauss]);

% Compute tensor products
rirj = drExp .* drTExp;
rinj = drExp .* permute(normExp, [2 1 3 4]);
rjni = normExp .* drTExp;

% Identity matrix
I = repmat(eye(D), [1, 1, numCtrlPts, numGauss]);

% Compute the tensor terms
Ttemp = 1 ./ (rExp.^2) * const3 .* (-drdnExp .* (const4 * I + D * rirj) + ...
    const4 * (rinj - rjni));
Utemp = const1 ./ rExp .* (rirj + I * const2);

% Gaussian weights
gwtUExp = repmat(reshape(gwtU(:), [1, 1, 1, numGauss]), [1, 1, numCtrlPts, 1]);
gwtVExp = repmat(reshape(gwtV(:), [1, 1, 1, numGauss]), [1, 1, numCtrlPts, 1]);

% Compute H and G matrices
H = jacobExp .* Ttemp .* gwtUExp .* gwtVExp .* NExp;
G = jacobExp .* Utemp .* gwtUExp .* gwtVExp .* NExp;

H = sum(H, 4);
G = sum(G, 4);

H = reshape(H, D, []);
G = reshape(G, D, []);
end
