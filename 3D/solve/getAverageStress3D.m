function avgStress = getAverageStress3D(nurbsStr, ctrlTrac)
% Calculate the average stress by integral along the external boundary.
%
% Inputs:
%   nurbsStr: Structure array containing boundary element data.
%   ctrlTrac: Traction values at control points.
%
% Outputs:
%   avgStress: Average stress vector.

avgStress = zeros(1, 6);
ngp_r = 16; % Number of Gauss points
[gpt, gwt] = lgwt(ngp_r, -1, 1);
[gptU, gptV] = meshgrid(gpt, gpt);
[gwtU, gwtV] = meshgrid(gwt, gwt);

D = 3;
indT = zeros(7, 1);
for i = 1:6
    indT(i + 1) = indT(i) + nurbsStr(i).numCollocT;
    for e = 1:nurbsStr(i).numElem
        range = nurbsStr(i).elemRange(e, :);
        gloIdx = nurbsStr(i).globIdx(e, :);

        % The coordinate of control points in the element
        elemCoords = nurbsStr(i).ctrlPts(gloIdx, 1:D);
        p = nurbsStr(i).p;
        q = nurbsStr(i).q;
        weights = nurbsStr(i).weights;
        tConn = nurbsStr(i).tConn(e, :);
        col1 = zeros(D * (p + 1) * (q + 1), 1);
        for k = 1:D
            col1(k:D:end) = tConn * D - D + k;
        end
        col1 = col1 + indT(i) * D;
        Ti = ctrlTrac(col1);
        Ti = [Ti(1:D:end), Ti(2:D:end), Ti(3:D:end)];

        uConn = nurbsStr(i).uConn(e, :);
        Xj = nurbsStr(i).ctrlPts(uConn, :);

        avgStressSub = integralAverageStress(elemCoords, gloIdx, range, p, q, ...
            nurbsStr(i).knotU, nurbsStr(i).knotV, weights, gptU, gptV, gwtU, gwtV, Ti, Xj);
        avgStress = avgStress + avgStressSub;
    end
end
end

function avgStressSub = integralAverageStress(elemCoords, gloIdx, range, p, q, ...
    knotU, knotV, weights, gptU, gptV, gwtU, gwtV, Ti, Xj)
% Computes the average stress contribution from a single element.

% Jacobian from parent to parameter space
jacobParam = (range(2) - range(1)) * (range(4) - range(3)) / 4;

% Gauss points in parameter space
xiU = convertToParamSpace(gptU(:), range(1:2));
xiV = convertToParamSpace(gptV(:), range(3:4));
numGauss = length(xiU);

% Compute NURBS basis functions and derivatives at Gauss points
N = zeros(1, (p + 1) * (q + 1), numGauss);
dNu = N;
dNv = N;
for i = 1:numGauss
    [N(:, :, i), dNu(:, :, i), dNv(:, :, i)] = NURBS2DBasisDers([xiU(i, :), xiV(i, :)], ...
        p, q, knotU, knotV, weights');
end

% Compute geometry derivatives
dxi = pagemtimes(dNu, elemCoords);
deta = pagemtimes(dNv, elemCoords);

% Compute normals
normals = cross(dxi, deta, 2);
jacobXi = vecnorm(normals, 2, 2);
jacob = jacobXi * jacobParam; % The final Jacobian we use

avgStressSub = computeStressKernel(numGauss, jacob, N, gwtU, gwtV, Ti, Xj);
end

function avgStressSub = computeStressKernel(numGauss, jacob, N, gwtU, gwtV, Ti, Xj)
% Calculate the average stress.
n = 9;

% Compute Ti * Xj for each stress component
TiXj11 = Ti(:, 1) * Xj(:, 1)';
TiXj22 = Ti(:, 2) * Xj(:, 2)';
TiXj33 = Ti(:, 3) * Xj(:, 3)';
TiXj12 = Ti(:, 1) * Xj(:, 2)';
TiXj21 = Ti(:, 2) * Xj(:, 1)';
TiXj23 = Ti(:, 2) * Xj(:, 3)';
TiXj32 = Ti(:, 3) * Xj(:, 2)';
TiXj13 = Ti(:, 1) * Xj(:, 3)';
TiXj31 = Ti(:, 3) * Xj(:, 1)';

% Symmetric parts
TiXj12M = (TiXj12 + TiXj21) / 2;
TiXj23M = (TiXj23 + TiXj32) / 2;
TiXj13M = (TiXj13 + TiXj31) / 2;

% Expand for Gauss points
TiXj11Exp = reshape(repmat(TiXj11, [1, 1, numGauss]), [1, n * n, numGauss]);
TiXj22Exp = reshape(repmat(TiXj22, [1, 1, numGauss]), [1, n * n, numGauss]);
TiXj33Exp = reshape(repmat(TiXj33, [1, 1, numGauss]), [1, n * n, numGauss]);
TiXj12Exp = reshape(repmat(TiXj12M, [1, 1, numGauss]), [1, n * n, numGauss]);
TiXj23Exp = reshape(repmat(TiXj23M, [1, 1, numGauss]), [1, n * n, numGauss]);
TiXj13Exp = reshape(repmat(TiXj13M, [1, 1, numGauss]), [1, n * n, numGauss]);

% NURBS basis functions expansion
NExp = reshape(N, [1, n, numGauss]);
Nij = pagemtimes(NExp, 'transpose', NExp, 'none');
NijExp = reshape(Nij, [1, n * n, numGauss]);

% Compute the Jacobian
jacob = reshape(jacob, [1, 1, numGauss]);

% Gaussian weights
gwtUExp = reshape(gwtU(:), [1, 1, numGauss]);
gwtVExp = reshape(gwtV(:), [1, 1, numGauss]);

% Compute stress components
H11 = jacob .* gwtUExp .* gwtVExp .* NijExp .* TiXj11Exp;
H22 = jacob .* gwtUExp .* gwtVExp .* NijExp .* TiXj22Exp;
H33 = jacob .* gwtUExp .* gwtVExp .* NijExp .* TiXj33Exp;
H12 = jacob .* gwtUExp .* gwtVExp .* NijExp .* TiXj12Exp;
H23 = jacob .* gwtUExp .* gwtVExp .* NijExp .* TiXj23Exp;
H13 = jacob .* gwtUExp .* gwtVExp .* NijExp .* TiXj13Exp;

H11 = sum(H11, 'all');
H22 = sum(H22, 'all');
H33 = sum(H33, 'all');
H12 = sum(H12, 'all');
H23 = sum(H23, 'all');
H13 = sum(H13, 'all');

avgStressSub = [H11, H22, H33, H12, H23, H13];
end
