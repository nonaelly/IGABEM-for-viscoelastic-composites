function avgStressSub = integralAvgStress2D(elemCoords, gloIdx, range, p, ...
        knotU, weights, gptU, gwtU, Ti, Xj)
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
    [N(:,:,index), dN(:,:,index)] = NURBSbasis(gloIdx, p, ...
        xiParamSpace(index,:), knotU, weights');
end

dxydxi = pagemtimes(dN, elemCoords); 
jacobXi = vecnorm(dxydxi, 2, 2); % Norm in every row
jacob = jacobXi * jacobParam; % The final jacobian we use

avgStressSub = AvgStress(numGauss, jacob, N, gwtU, Ti, Xj);

end

%% Sub functions.
% =========================================================================
% Calculate the average stress.
function avgStress = AvgStress(numGauss, jacob, N, gwtU, Ti, Xj)

n = 3;

TiXj11 = Ti(:, 1) *Xj(:, 1)';
TiXj22 = Ti(:, 2) *Xj(:, 2)';
TiXj12 = Ti(:, 1) *Xj(:, 2)';
TiXj21 = Ti(:, 2) *Xj(:, 1)';
TiXj12M = (TiXj12 + TiXj21)/2;

TiXj11Exp = repmat(TiXj11, [1, 1, numGauss]);
TiXj11Exp = reshape(TiXj11Exp, [1, n*n, numGauss]);

TiXj22Exp = repmat(TiXj22, [1, 1, numGauss]);
TiXj22Exp = reshape(TiXj22Exp, [1, n*n, numGauss]);

TiXj12Exp = repmat(TiXj12M, [1, 1, numGauss]);
TiXj12Exp = reshape(TiXj12Exp, [1, n*n, numGauss]);

% n = 3, represent 3 control points.
NExp = reshape(N, [1, n, numGauss]);

% Nij: 3*3
Nij = pagemtimes(NExp,'transpose', NExp, 'none');
NijExp = reshape(Nij, [1, n*n, numGauss]);

jacob = reshape(jacob, [1, 1, numGauss]);

arrayGwtU = reshape(gwtU(:), [1, 1, numGauss]);

NijTiXj11 = NijExp .* TiXj11Exp;
NijTiXj22 = NijExp .* TiXj22Exp;
NijTiXj12 = NijExp .* TiXj12Exp;

S11 = jacob .*arrayGwtU .*NijTiXj11;
S22 = jacob .*arrayGwtU .*NijTiXj22;
S12 = jacob .*arrayGwtU .*NijTiXj12;

S11 = sum(S11, 'all');
S22 = sum(S22, 'all');
S12 = sum(S12, 'all');

avgStress = [S11 S22 S12];

end