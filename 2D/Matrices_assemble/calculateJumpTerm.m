function [Cij] = calculateJumpTerm(n1, n2, nu)
% Calculate the jump term at a collocation point using the formula from the paper
% Inputs:
%   n1: Normal vector of the first element
%   n2: Normal vector of the second element
%   mu: Shear modulus
% Output:
%   Cij: Jump term matrix

% Initialize
Cij = zeros(2, 2);
t = zeros(2, 2);

% Calculate tangential vectors
tempTang = cross([n1 0], [0 0 1]);
t(1, :) = tempTang(1:2);
tempTang = cross([n2 0], [0 0 1]);
t(2, :) = tempTang(1:2);

% Calculate dot product and angles
dotProd = dot(n1, n2);
if (dotProd - 1) >= eps
    theta = 0;
else
    theta = acos(dot(n1, n2));
end

if n1(1) * n2(2) > n1(2) * n2(1)
    thetaBar = pi - theta;
else
    thetaBar = pi + theta;
end

xVector = [1 0];
elmTheta = zeros(1, 2);

for i = 1:2
    theta = acos(dot(xVector, t(i, :)));
    if xVector(1) * t(i, 2) > xVector(2) * t(i, 1)
        elmTheta(i) = theta;
    else
        elmTheta(i) = 2 * pi - theta;
    end
end

% Calculate the jump term using the formula
term1 = 1 / (8 * pi * (1 - nu));
term2 = 4 * (1 - nu) * thetaBar;
term3 = sin(2 * elmTheta(1)) - sin(2 * elmTheta(2));
term4 = cos(2 * elmTheta(2)) - cos(2 * elmTheta(1));

Cij(1, 1) = (term2 + term3) * term1;
Cij(1, 2) = term1 * term4;
Cij(2, 1) = term1 * term4;
Cij(2, 2) = (term2 - term3) * term1;

end
