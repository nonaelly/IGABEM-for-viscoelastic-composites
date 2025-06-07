function [T, U] = getSingularityPart(xiQ, xiP, elCoords, collocCoords, p, q, ...
    knotU, knotV, weights, const1, const2, const3, const4, id)
% Compute the singularity part for strong singular integrals using the power series method
%
% Inputs:
%   xiQ, xiP: Parametric coordinates of the integration and collocation points
%   elCoords: Element coordinates
%   collocCoords: Coordinates of the collocation points
%   p, q: Degrees of the NURBS basis functions
%   knotU, knotV: Knot vectors in the U and V directions
%   weights: Weights of the NURBS basis functions
%   const1, const2, const3, const4: Material constants
%   id: Element ID
%
% Outputs:
%   T, U: Strong singular part of the T and U matrices

% Initial calculations
xiDiff = xiQ - xiP;
rhoQp = norm(xiDiff);
rhoQpXi = xiDiff(1) / rhoQp;
rhoQpEta = xiDiff(2) / rhoQp;

% Compute derivatives of the NURBS basis functions
[~, dNu, dNv] = NURBS2DBasisDers(xiP, p, q, knotU, knotV, weights);
dxi = dNu * elCoords;
deta = dNv * elCoords;
drho = dxi * rhoQpXi + deta * rhoQpEta;

M = 6; % Number of terms in the power series
G = zeros(M + 1, 1);
G(1) = norm(drho)^2; % G0

rhoM = linspace(0, rhoQp, M + 1);
RM = rhoM(2:end)'.^(0:M-1);
rhoMXi = linspace(xiP(1), xiQ(1), M + 1);
rhoMEta = linspace(xiP(2), xiQ(2), M + 1);
rM = zeros(1, M);

for row = 1:M
    [N, ~, ~] = NURBS2DBasisDers([rhoMXi(row+1), rhoMEta(row+1)], p, q, knotU, knotV, weights);
    rM(row) = norm(N * elCoords - collocCoords);
end

Z = ((rM ./ rhoM(2:end)).^2 - G(1)) ./ rhoM(2:end);
G(2:end) = RM \ Z';

% Coefficients C
C = zeros(M + 1, 1);
C(1) = sqrt(G(1)); % C0
for i = 1:M
    temp = 0;
    for j = 1:i-1
        temp = temp + C(j + 1) * C(i - j + 1);
    end
    C(i + 1) = 1 / (2 * C(1)) * (G(i + 1) - temp);
end

% Coefficients H
H = zeros(5, 1);
H(1) = 1 / C(1);
CBar = C / C(1);
H(2) = CBar(2);
H(3) = 2 * CBar(3) - CBar(2)^2;
H(4) = 3 * CBar(4) - 3 * CBar(2) * CBar(3) + CBar(2)^3;
H(5) = 4 * CBar(5) + 4 * CBar(2)^2 * CBar(3) - 4 * CBar(2) * CBar(4) - 2 * CBar(3)^2 - CBar(2)^4;

% Calculation fBar.
K = M;
rhoKXi = linspace(xiP(1), xiQ(1), K + 1);
rhoKEta = linspace(xiP(2), xiQ(2), K + 1);
jacobians = zeros(1, K);
Tbar = zeros(K + 1, 3, 3 * (p + 1) * (q + 1));
Ubar = zeros(K + 1, 3, 3 * (p + 1) * (q + 1));

for row = 1:K + 1
    [N, dNu, dNv] = NURBS2DBasisDers([rhoKXi(row), rhoKEta(row)], p, q, knotU, knotV, weights);
    dxdxi = dNu * elCoords;
    dxdeta = dNv * elCoords;
    n = cross(dxdxi, dxdeta, 2);
    jacobians(row) = vecnorm(n, 2, 2);

    % The norm of boundary inside.
    if id > 6
        n = -n ./ jacobians(row);
    else
        n = n ./ jacobians(row);
    end

    if row == 1
        dr = drho / norm(drho);
        drdn = dr * n';
    else
        fieldPt = pagemtimes(N, elCoords);
        relDist = fieldPt - collocCoords;
        rK = vecnorm(relDist, 2, 2);
        dr = 1 ./ rK .* relDist;
        drdn = dr * n';
    end

    % Compute singularity part
    T = zeros(3, 3);
    U = zeros(3, 3);
    for i = 1:3
        for j = 1:3
            [T(i, j), U(i, j)] = getKernelPSE(i, j, dr, drdn, n, const1, const2, const3, const4);
        end
    end

    for i = 1:(p + 1) * (q + 1)
        Tbar(row, :, i * 3 - 2:i * 3) = N(i) * T;
        Ubar(row, :, i * 3 - 2:i * 3) = N(i) * U;
    end
end

rhoN = linspace(0, rhoQp, K + 1);
RN = rhoN(2:end)'.^(0:K-1);
rhoBar = zeros(K, 1);
for i = 2:K + 1
    rhoBar(i - 1) = sqrt((rhoN(i).^(0:M)) * G);
end

lambda = 2;
T = getF(Tbar, K, lambda, rhoBar, jacobians, rhoN, RN, C, rhoQp, H, p, q);
lambda = 1;
U = getF(Ubar, K, lambda, rhoBar, jacobians, rhoN, RN, C, rhoQp, H, p, q);
end

function [Tij, Uij] = getKernelPSE(i, j, dr, drdn, normal, const1, const2, const3, const4)
% Compute the kernel for the power series expansion
%
% Inputs:
%   i, j: Indices for the kernel
%   dr, drdn: Derivatives of r and n
%   normal: Normal vector
%   const1, const2, const3, const4: Material constants
%
% Outputs:
%   Tij, Uij: Kernel values

Tij = const3 * (-drdn * (const4 * (i == j) + 3 * (dr(i) * dr(j))) ...
    + const4 * (dr(i) * normal(j) - dr(j) * normal(i)));

Uij = const1 * (dr(i) * dr(j) + (i == j) * const2);
end

function F = getF(fBar, N, lamda, rhoBar, jacobians, rhoN, RN, C, rhoqp, H, p, q)
% Compute the F matrix using the power series expansion
%
% Inputs:
%   fBar: Precomputed F-bar values
%   N: Number of terms in the series
%   lambda: Lambda parameter
%   rhoBar: Rho-bar values
%   jacobians: Jacobian values
%   rhoN: Rho values
%   RN: Rho N values
%   C: Coefficient C values
%   rhoQp: Rho Qp value
%   H: Coefficient H values
%   p, q: Degrees of the NURBS basis functions
%
% Outputs:
%   F: Computed F matrix

B = zeros(N + 1, 3, 3 * (p + 1) * (q + 1));
B(1, :, :) = fBar(1, :, :) * jacobians(1) / C(1)^lamda;

Y = zeros(N, 3, 3 * (p + 1) * (q + 1));
for i = 1:N
    Y(i, :, :) = (permute(fBar(i + 1, :, :), [2, 3, 1]) * jacobians(i + 1) / rhoBar(i).^lamda - permute(B(1, :, :), [2, 3, 1])) / rhoN(i + 1);
end
for i = 1:3
    for j = 1:3 * (p + 1) * (q + 1)
        B(2:end, i, j) = RN \ Y(:, i, j);
    end
end

F = 0;  % [3,3*(p+1)*(q+1)]
for n = 0:N
    if 0 <= n && n <= lamda-3
        G = 1/(n-lamda+2)*(1/rhoqp^(lamda-n-2)-H(lamda-n-2+1));
    elseif n == lamda-2
        G = log(rhoqp) - log(H(1));
    else
        G = rhoqp^(n-lamda+2)/(n-lamda+2);
    end
    F = F + permute(B(n+1, :, :), [2, 3, 1]) * G;
end
end
