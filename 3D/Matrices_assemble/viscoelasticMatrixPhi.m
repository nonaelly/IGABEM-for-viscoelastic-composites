function rimPhi = viscoelasticMatrixPhi(approxPts, appliedPts, dA, dim)
% Function: viscoelasticMatrixPhi
% Description: Forms the transformation matrix for the radial integration method (RIM).
%
% Input:
%   approxPts: Points to be approximated (N x D matrix).
%   appliedPts: Applied points in RIM, including boundary and inner points (M x D matrix).
%   dA: Size of the support domain.
%   dim: Dimension (2 or 3).
%
% Output:
%   rimPhi: Transformation matrix in radial integration method.
%
% Note: For example, u = rimPhi(1:end-2*(dim+1), :) * alpha

% -------------------------------------------------------------------------
numApplied = size(appliedPts, 1);
numApprox = size(approxPts, 1);
rimPhi = zeros((numApprox + dim + 1) * dim, (numApplied + dim + 1) * dim);

for i = 1:numApprox
    rowIdx = i * dim - dim + 1 : i * dim;
    for j = 1:numApplied
        colIdx = j * dim - dim + 1 : j * dim;
        R = norm(approxPts(i,:) - appliedPts(j,:));
        if R < dA
            phiR = 1 - 6 * ((R / dA)^2) + 8 * ((R / dA)^3) - 3 * ((R / dA)^4);
        else
            phiR = 0;
        end
        phiSubMatrix = phiR * eye(dim);
        rimPhi(rowIdx, colIdx) = phiSubMatrix;
    end
end

% Handling the identity part
colIdx = (numApplied + 1) * dim - dim + 1 : (numApplied + 1) * dim;
for i = 1:numApprox
    rowIdx = i * dim - dim + 1 : i * dim;
    phiSubMatrix = eye(dim);
    rimPhi(rowIdx, colIdx) = phiSubMatrix;
end

% Handling the linear part
for i = 1:numApprox    
    rowIdx = i * dim - dim + 1 : i * dim;
    for j = numApplied + 2 : numApplied + dim + 1
        colIdx = j * dim - dim + 1 : j * dim;
        xk = approxPts(i, j - 1 - numApplied);
        phiSubMatrix = xk * eye(dim);
        rimPhi(rowIdx, colIdx) = phiSubMatrix;
    end
end

% Handling the applied points part
rowIdx = (numApprox * dim) + 1 : (numApprox + dim + 1) * dim;
for j = 1:numApplied
    colIdx = j * dim - dim + 1 : j * dim;
    x1 = appliedPts(j, 1);
    x2 = appliedPts(j, 2);
    phiSubMatrix = [eye(dim); x1 * eye(dim); x2 * eye(dim)];
    if dim == 3
        x3 = appliedPts(j, 3);
        phiSubMatrix = [phiSubMatrix; x3 * eye(dim)];
    end
    rimPhi(rowIdx, colIdx) = phiSubMatrix;
end

end
