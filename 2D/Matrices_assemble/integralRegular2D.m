function [H, G] = integralRegular2D(elemCoords, gloBasisIdx, collocCoords, ...
    range, p, knotVec, weights, const1, const2, const3, const4, gptU, gwtU)
% integralTelles2D: Computes the G matrix for the weakly singular case using
% the Telles transformation.
%
% Inputs:
%   - elemCoords: Coordinates of the element.
%   - gloBasisIdx: Global basis functions with non-zero support in the element.
%   - collocCoords: Coordinates of the collocation points.
%   - xiParam: Parameter coordinate of the collocation point.
%   - range: Range of the element in parameter space [xi_i, xi_i+1].
%   - p: Degree of the NURBS curve.
%   - knotVec: Knot vector of the NURBS curve.
%   - weights: Weights of the control points.
%   - const1, const2, const3, const4: Constants used in the calculation.
%   - gptU: Gauss points in the parameter space.
%   - gwtU: Gauss weights in the parameter space.
%
% Output:
%   - G: The G matrix for the weakly singular case.

% -------------------------------------------------------------------------
% Jacobian from parent to parameter space
jacobParam = (range(2) - range(1)) / 2;

% Gauss points in parameter space
xiParamSpace = convertToParamSpace(gptU(:), range);

numGauss = length(xiParamSpace); % Number of Gauss points

% Preallocate arrays
N = zeros(1, p + 1, numGauss); % 'N_Gauss' pages. For each page, 1 * (p + 1)
dN = N;

% Compute NURBS basis functions and derivatives at Gauss points
for index = 1:numGauss
    [N(:,:,index), dN(:,:,index)] = NURBSbasis(gloBasisIdx, p, ...
        xiParamSpace(index,:), knotVec, weights');
end

[jacobXi, normals, r, dr, drdn] = getKernelParameters(elemCoords, ...
    collocCoords, N, dN);

jacob = jacobXi * jacobParam; % The final jacobian we use

% Compute the G matrix using ElastKernel2D
[H, G] = ElastKernel2D(r, dr, drdn, const1, const2, const3, const4, numGauss, ...
    jacob, normals, N, gwtU, 2);

end

%% Sub functions.
% =========================================================================
% Calculate the regular matrix.
function [H, G] = ElastKernel2D(r, dr, drdn, const1, const2, const3, const4, ...
    numGauss, jacob, normals, N, gwtU, D)
% ElastKernelTelles2D: Calculates the elastic kernel matrix using the
% Telles transformation for weakly singular integrals.
%
% Inputs:
%   - r: Distance between field points and source coordinates.
%   - dr: Derivative of r.
%   - drdn: Dot product of dr and normals.
%   - const1, const2, const3, const4: Constants for elasticity calculation.
%   - numGauss: Number of Gaussian points.
%   - jacob: Jacobian of the transformation.
%   - normals: Normal vectors.
%   - N: NURBS basis functions.
%   - gwtU: Gaussian weights.
%   - D: Dimensionality of the problem.
%
% Outputs:
%   - G: Elastic kernel matrix.

numCtrlPts = size(N, 2); % Number of control points

% Reshape matrices for element-wise operations
r = repmat(r, [numCtrlPts, 1, 1]);
rExp = reshape(r, [1, 1, numCtrlPts, numGauss]);

dr = repmat(dr, [numCtrlPts, 1, 1]);
dr = permute(dr, [2, 1, 3]);
drExp = reshape(dr, [D, 1, numCtrlPts, numGauss]);
drTExp = permute(drExp, [2, 1, 3, 4]);

drdn = repmat(drdn, [numCtrlPts, 1, 1]);
drdnExp = reshape(drdn, [1, 1, numCtrlPts, numGauss]);

normals = repmat(normals, [numCtrlPts, 1, 1]);
normals = permute(normals, [2, 1, 3]);
normExp = reshape(normals, [D, 1, numCtrlPts, numGauss]);
normTExp = permute(normExp, [2, 1, 3, 4]);

NExp = reshape(N, [1, 1, numCtrlPts, numGauss]);

% Compute tensor products
riRj = drExp .* drTExp;
riNj = drExp .* normTExp;
rjNi = normExp .* drTExp;

% Identity matrix
I0 = eye(D);
I = repmat(I0, [1, 1, numCtrlPts, numGauss]);

% Reshape and repeat matrices for element-wise operations
jacob = reshape(jacob, [1, 1, 1, numGauss]);
jacobExp = repmat(jacob, [1, 1, numCtrlPts, 1]);

gwtUExp = reshape(gwtU(:), [1, 1, 1, numGauss]);
gwtExp = repmat(gwtUExp, [1, 1, numCtrlPts, 1]);

% Compute the tensor term
Ttemp = 1 ./ rExp .* const3 .* (-drdnExp .* (const4 * I + 2 * riRj) + const4 * (riNj - rjNi));
Utemp = const1 * (riRj + I .* (const2 * log(1 ./ rExp)));

% Compute the elastic kernel matrix
H2 = jacobExp .* gwtExp .* Ttemp .* NExp;
G2 = jacobExp .* gwtExp .* Utemp .* NExp;

H2 = sum(H2, 4);
G2 = sum(G2, 4);

H = reshape(H2, D, []);
G = reshape(G2, D, []);

end
