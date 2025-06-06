function H = integralSST2D(elemCoords, gloBasisIdx, collocCoords, xiParam, ...
    range, p, knotVec, weights, Cij, ~, ~, const3, const4, gptU, gwtU)
% integralSST2D: Computes the H matrix for the singular case using the
%               subtraction of singularity technique by Guiggiani and Casalini.
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
%   - Cij: Jump term for singularities.
%   - const1, const2, const3, const4: Constants used in the calculation.
%   - gptU: Gauss points in the parameter space.
%   - gwtU: Gauss weights in the parameter space.
%   - D: Dimension of the problem (2 for 2D).
%   - id: Element ID.
%
% Output:
%   - H: The H matrix for the singular case.

% -------------------------------------------------------------------------
D = 2;

% Fix the parameter coordinate to avoid singularities
if xiParam == range(1)
    nudgedParam = xiParam + eps;
elseif xiParam == range(2)
    nudgedParam = xiParam - eps;
else
    nudgedParam = xiParam;
end

% Get the coordinate in the Gauss space [-1, 1]
collocCoord = convertToParentCoordSpace(xiParam, range);

% Calculate NURBS basis functions and derivatives at the nudged parameter
[srcN, srcdN] = NURBSbasis(gloBasisIdx, p, nudgedParam, knotVec, weights');

% Construct htermMatrix
hterm = [0, -const4 * const3; const4 * const3, 0];
htermMatrix = [hterm .* srcN(1), hterm .* srcN(2), hterm .* srcN(3)];

% Compute Jacobian from parent to parameter space
jacobParam = (range(2) - range(1)) / 2;

% Gauss points in parameter space
xiParamSpace = convertToParamSpace(gptU(:), range);

numGauss = length(xiParamSpace); % Number of Gauss points

% Preallocate arrays
N = zeros(1, p+1, numGauss); % 'N_Gauss' pages, each page 1*p+1
dN = N;

% Compute NURBS basis functions and derivatives at Gauss points
for i = 1:numGauss
    [N(:, :, i), dN(:, :, i)] = NURBSbasis(gloBasisIdx, p, xiParamSpace(i, :), knotVec, weights');
end

% Get kernel parameters
[jacobXi, normals, r, dr, drdn] = getKernelParameters(elemCoords, collocCoords, N, dN);

% Compute final Jacobian
jacob = jacobXi * jacobParam;

% Compute det_xi
detXi = gptU(:) - collocCoord;
detXi = reshape(detXi, [1, 1, numGauss]);

% Compute H matrix using ElastKernelSST2D
H = ElastKernelSST2D(r, dr, drdn, const3, const4, numGauss, jacob, normals, N, gwtU, D, detXi, htermMatrix, p);

% Add analytical integral
if abs(collocCoord) < (1 - 100 * eps)
    % collocCoord belongs to (-1, 1)
    H = H + htermMatrix * log(abs((1 - collocCoord) / (1 + collocCoord)));
else
    % collocCoord = -1/1
    jacobS = getKernelParameters(elemCoords, collocCoords, srcN, srcdN) * jacobParam;
    H = H + htermMatrix * log(abs(jacobS)) * sign(-collocCoord);
end

% Add jump terms
if abs(xiParam - range(2)) > 100 * eps
    jumpMatrix = [Cij .* srcN(1), Cij .* srcN(2), Cij .* srcN(3)];
    H = H + jumpMatrix;
end

end

%% Sub functions.
% =========================================================================
% Calculate the strong singular matrix.
function H = ElastKernelSST2D(r, dr, drdn, const3, const4, numGauss, jacob, ...
    normals, N, gwtU, D, detXi, htermMatrix, p)
% ElastKernelSST2D: Calculates the elastic kernel matrix using the
%   subtraction of singularity technique (SST).
%
% Inputs:
%   - r: Distance between field points and source coordinates.
%   - dr: Derivative of r.
%   - drdn: Dot product of dr and normals.
%   - const3, const4: Constants for elasticity calculation.
%   - numGauss: Number of Gaussian points.
%   - jacob: Jacobian of the transformation.
%   - normals: Normal vectors.
%   - N: NURBS basis functions.
%   - gwtU: Gaussian weights.
%   - D: Dimensionality of the problem.
%   - detXi: xi-xi'.
%   - htermMatrix: Matrix for handling singularity.
%
% Outputs:
%   - H: Elastic kernel matrix.

numCtrlPts = p+1; % Number of control points

% Reshape matrices for element-wise operations
r = repmat(r, [p+1, 1, 1]);
rExp = reshape(r, [1, 1, numCtrlPts, numGauss]);

dr = repmat(dr, [p+1, 1, 1]);
dr = permute(dr, [2, 1, 3]);
drExp = reshape(dr, [D, 1, numCtrlPts, numGauss]);
drTExp = permute(drExp, [2, 1, 3, 4]);

drdn = repmat(drdn, [p+1, 1, 1]);
drdnExp = reshape(drdn, [1, 1, numCtrlPts, numGauss]);

normals = repmat(normals, [p+1, 1, 1]);
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
detXi = reshape(detXi, [1, 1, 1, numGauss]);

gwtUExp = reshape(gwtU(:), [1, 1, 1, numGauss]);
gwtExp = repmat(gwtUExp, [1, 1, numCtrlPts, 1]);

% Compute the tensor term
Ttemp = 1 ./ rExp .* const3 .* (-drdnExp .* (const4 * I + 2 * riRj) + const4 * (riNj - rjNi));

% Reshape and repeat matrices for element-wise operations
hExp = repmat(htermMatrix, [1, 1, numGauss]);
hExp = reshape(hExp, [D, D, numCtrlPts, numGauss]);

% Compute the elastic kernel matrix
H1 = (jacobExp .* gwtExp .* Ttemp .* NExp .* detXi - hExp) ./ detXi;
H1 = sum(H1, 4);
H = reshape(H1, D, []);

end
