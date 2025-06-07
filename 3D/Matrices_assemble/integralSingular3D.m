function [H, G] = integralSingular3D(ngp, xi, eta, p, q, range, elCoords, ...
    collocCoords, knotU, knotV, weights, const1, const2, const3, const4, id)
% integralSingular3D: Computes the strong singular integrals using the power series method.
%
% Inputs:
%   ngp: Number of Gaussian points.
%   xi_u, xi_v: Parametric coordinates of the collocation point.
%   p, q: Degrees of the NURBS basis functions.
%   flag: Array indicating whether the collocation point is on each element side.
%   range: Parametric range of the element.
%   elCoords: Coordinates of the element control points.
%   collocCoords: Coordinates of the collocation points.
%   knotVecU, knotVecV: Knot vectors in the U and V directions.
%   weights: Weights of the NURBS basis functions.
%   const1, const2, const3, const4: Material constants.
%   id: Element ID.
%
% Outputs:
%   H: Strong singular part of the H matrix.
%   G: Strong singular part of the G matrix.
% 
%         3             
%      D-----C     
%      |     |      
%    4 |     | 2   
%      |     |       
%      A-----B       
%         1   

% Gauss points and weights
[gpt, gwt] = lgwt(ngp, -1, 1);
xip = [xi, eta];
H = zeros(3, 3 * (p + 1) * (q + 1));
G = zeros(3, 3 * (p + 1) * (q + 1));

% u 4/2 correspond -1/1
% v 1/3 correspond -1/1
flag = getFlagPSE(xi, eta, range);

for sideIndex = 1:4
    if flag(sideIndex)
        continue;
    end
    
    switch sideIndex
        case 1
            gaussDir = 1; % High Gaussian integration in the 'u' direction
            rangeIdx = 3;
            normalVec = [0, -1];
            tol = 1e-15;
        case 2
            gaussDir = 2; % High Gaussian integration in the 'v' direction
            rangeIdx = 2;
            normalVec = [1, 0];
            tol = -1e-15;
        case 3
            gaussDir = 1; % High Gaussian integration in the 'u' direction
            rangeIdx = 4;
            normalVec = [0, 1];
            tol = -1e-15;
        case 4
            gaussDir = 2; % High Gaussian integration in the 'v' direction
            rangeIdx = 1;
            normalVec = [-1, 0];
            tol = 1e-15;
    end
    
    contDir = 3 - gaussDir; % Continuous direction perpendicular to Gaussian integration direction
    Jacob = (range(gaussDir * 2) - range(gaussDir * 2 - 1)) / 2;
    
    for pt = 1:ngp
        xiq = zeros(1, 2);
        xiq(contDir) = range(rangeIdx) + tol;
        xiq(gaussDir) = convertToParamSpace(gpt(pt), range(gaussDir * 2 - 1:gaussDir * 2));
        xiqp = xiq - xip;
        rhoqp = norm(xiqp);
        drhodn = xiqp * normalVec' / rhoqp;
        
        [FH, FG] = getSingularityPart(xiq, xip, elCoords, collocCoords, ...
            p, q, knotU, knotV, weights, const1, const2, const3, const4, id);
        
        H = H + FH / rhoqp * drhodn * Jacob * gwt(pt);
        G = G + FG / rhoqp * drhodn * Jacob * gwt(pt);
    end
end
end
