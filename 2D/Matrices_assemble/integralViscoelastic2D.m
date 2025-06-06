function [HVis1_1, HVis2_1, HVis1_2, HVis2_2] = integralViscoelastic2D(elemCoords, ...
    gloBasisIdx, sourcePt, xi, range, p, knotVec, weights, const1, const2, const3, ...
    const4, gptU_r, gptU_s, gwtU_r, gwtU_s, appliedPts, dA, flag, D)
% integralViscoelastic2D: Computes the viscoelastic matrix H for both
%   weak singular and regular cases using the analysis solution in RIM.
%
% Inputs:
%   - elemCoords: Coordinates of the element.
%   - gloBasisIdx: Global basis functions.
%   - sourcePt: Coordinates of the source points.
%   - xi: Source parameter coordinate.
%   - range: Range of the element in parameter space.
%   - p: Degree of the NURBS curve.
%   - knotVec: Knot vector of the NURBS curve.
%   - weights: Weights of the control points.
%   - const1, const2, const3, const4: Constants for elasticity calculation.
%   - gptU_r: Gauss points in the parameter space (Regular).
%   - gwtU_r: Gauss weights in the parameter space (Regular).
%   - gptU_s: Gauss points in the parameter space (Singular).
%   - gwtU_s: Gauss weights in the parameter space (Singular).
%   - appliedPts: Applied points in RIM, containing boundary points and inner points.
%   - dA: Support domain size in RIM.
%   - flag: Flag for singularity.
%   - D: Dimensionality of the problem.
%
% Outputs:
%   - HVis1_1, HVis2_1: Viscoelastic matrix (column: 1 ~ num_applied).
%   - HVis1_2, HVis2_2: Viscoelastic matrix (column: 1+num_applied ~ 1+D+num_applied).

numApplied = size(appliedPts, 1);

% Jacobian from parent to parameter space
jacobParam = (range(2) - range(1)) / 2;

if flag
    srcXi = convertToParentCoordSpace(xi, range);
    xiStar = srcXi^2 - 1;
    gamBar = nthroot(srcXi * xiStar + abs(xiStar), 3) + ...
        nthroot(srcXi * xiStar - abs(xiStar) + srcXi, 3);
    xi = ((gptU_s - gamBar).^3 + gamBar * (gamBar.^2 + 3)) / (1 + 3 * gamBar.^2);
    xiParam = convertToParamSpace(xi, range);
    gwt = gwtU_s;
    jacobTelles = 3 * ((gptU_s - gamBar).^2) / (1 + 3 * gamBar.^2);
    jacobTelles = reshape(jacobTelles', [1, 1, size(gptU_s, 1)]);
else
    xiParam = convertToParamSpace(gptU_r(:), range);
    gwt = gwtU_r;
end

numGauss = length(xiParam); % Number of Gauss points (X and Y)

N = zeros(1, p + 1, numGauss); % 'N_Gauss' pages. For each page, 1*(p+1)
dN = N;

for index = 1:numGauss
    [N(:,:,index), dN(:,:,index)] = NURBSbasis(gloBasisIdx, p, ...
        xiParam(index,:), knotVec, weights');
end

[jacobXi, normals, r, dr, drdn] = getKernelParameters(elemCoords, sourcePt, N, dN);

jacob1 = jacobXi * jacobParam;

% The final jacobian we use
if flag
    jacob = jacob1 .* jacobTelles;
else
    jacob = jacob1;
end

[HVis1_1, HVis2_1] = ViscoKernelA2D(r, dr, drdn, const1, const4, appliedPts, ...
    sourcePt, numApplied, dA, numGauss, jacob, gwt, D);

%-----------------------------------------------------------
% The term without singularity.
gwtB = gwtU_r;
xiParamB = convertToParamSpace(gptU_r(:), range);

numGaussB = length(xiParamB); % Number of Gauss points

NB = zeros(1, p + 1, numGaussB); % 'N_Gauss' pages. For each page, 1*(p+1)
dNB = NB;

for index = 1:numGaussB
    [NB(:,:,index), dNB(:,:,index)] = NURBSbasis(gloBasisIdx, p, ...
        xiParamB(index,:), knotVec, weights');
end

[jacobXiB, ~, rB, drB, drdnB] = getKernelParameters(elemCoords, sourcePt, NB, dNB);

jacobB = jacobXiB * jacobParam; % The final jacobian we use

[HVis1_2, HVis2_2] = ViscoKernelB2D(rB, drB, drdnB, const1, const4, numGaussB, jacobB, gwtB, D);

end

%% Sub functions.
% =========================================================================
% HVis1_1, HVis2_1.
function [HVis1, HVis2] = ViscoKernelA2D(r, dr, drdn, const1, const4, appliedPts, ...
    sourcePt, numApplied, dA, numGauss, jacob, gwt, D)
% ViscoKernelA2D: Calculates the viscoelastic kernel matrices HVis1 and HVis2
% for weakly singular integrals using the analytical solution in RIM.
%
% Inputs:
%   - r, dr, drdn: Kernel parameters.
%   - const1, const4: Elasticity constants.
%   - appliedPts: Applied points in RIM, including boundary and inner points.
%   - sourcePt: Source collocation point.
%   - numApplied: Number of applied points.
%   - dA: Support domain size in RIM.
%   - numGauss: Number of Gauss points.
%   - jacob: Jacobian matrix.
%   - gwt: Gauss weights.
%   - D: Dimensionality of the problem (2 for 2D).
%
% Outputs:
%   - HVis1, HVis2: Viscoelastic matrices H1 and H2.

r = repmat(r, [numApplied, 1, 1]);
dr = repmat(dr, [numApplied, 1, 1]);
drdn = repmat(drdn, [numApplied, 1, 1]);
dr = permute(dr, [2, 1, 3]);

% Vectorized operations
collocGlbPtExp = repmat(sourcePt', [1, numApplied, numGauss]);
RBar = collocGlbPtExp - appliedPts';
RBarNorm = vecnorm(RBar, 2, 1);  % Compute the norm of each row vector
s = sum(dr .* RBar, 1);
RadialInte = zeros(D, numApplied * numGauss);

% Calculate the radial integral using analytical solution
% Case 1: RBarNorm.^2 - s.^2 > dA^2, radial integral = 0.

% Case 2: RBarNorm.^2 - s.^2 < dA^2
Jud2 = RBarNorm.^2 - s.^2 < dA^2;
GloInd2 = find(Jud2);
if any(Jud2(:))
    j1 = sqrt(dA^2 - RBarNorm(Jud2).^2 + s(Jud2).^2);
    r1 = -s(Jud2) - j1;
    r2 = -s(Jud2) + j1;
    % [a,b]
    a = max(0, r1);
    b = min(r2, r(Jud2));

    % Case 2.1: [a,b] is a empty assemble.

    % Case 2.2: [a,b] is a not empty assemble.
    Jud2_2 = a < b;
    LocInd2_2 = find(Jud2_2);
    GloInd2_2 = GloInd2(Jud2_2);
    if any(Jud2_2)
        % Case 2.2.1: RBar + s = 0 (A, p and Q are all in a line).
        Jud2_2_1 = abs(RBarNorm(GloInd2_2) + s(GloInd2_2)) < 1e-10;
        
        if any(Jud2_2_1)
            GloInd2_2_1 = GloInd2_2(Jud2_2_1);
            LocInd2_2_1 = find(Jud2_2_1);
            RadialInte(:, GloInd2_2_1) = computeCollinear(RBarNorm(:,GloInd2_2_1), ...
                s(:,GloInd2_2_1), dA, dr(:,GloInd2_2_1), RBar(:, GloInd2_2_1), D, ...
                a(LocInd2_2(LocInd2_2_1)), b(LocInd2_2(LocInd2_2_1)));
        end

        % Case 2.2.2: A, p and Q are not in a line.
        if any(~Jud2_2_1)
            GloInd2_2_2 = GloInd2_2(~Jud2_2_1);
            LocInd2_2_2 = find(~Jud2_2_1);
            RadialInte(:, GloInd2_2_2) = computeNonCollinear(RBarNorm(:,GloInd2_2_2), ...
                s(:,GloInd2_2_2), dA, dr(:,GloInd2_2_2), RBar(:, GloInd2_2_2), D, ...
                a(LocInd2_2(LocInd2_2_2)), b(LocInd2_2(LocInd2_2_2)));
        end
    end
end

RTempExp = reshape(RadialInte, D, numApplied, numGauss);
RExp = reshape(RTempExp, [D, 1, numApplied, numGauss]);
drExp = reshape(dr,[D, 1, numApplied, numGauss]);
drTExp = permute(drExp, [2 1 3 4]);

UiRj = RExp.*drTExp;
UjRi = permute(UiRj, [2 1 3 4]);
RiRj = drExp.*drTExp;
RkUk = sum(RExp.*drExp, 1);
I0 = eye(D);
I = repmat(I0, [1, 1, numApplied, numGauss]);

coeffient = drdn./(r.^(D-1));
coeffient = permute(coeffient, [2 1 3]);
coeffientExp = reshape(coeffient, [1, 1, numApplied, numGauss]);

jacob = reshape(jacob, [1, 1, 1, numGauss]);
jacobExp = repmat(jacob, [1, 1, numApplied, 1]);

Ttemp1 = -const1*coeffientExp.*( const4*(UiRj+RkUk.* I) - UjRi ...
    + D *RkUk.*RiRj );
Ttemp2 = -2*const4*const1*coeffientExp.* UjRi;

gwt = reshape(gwt(:), [1, 1, 1, numGauss]);
gwtExp = repmat(gwt, [1, 1, numApplied, 1]);

H1 = jacobExp.*Ttemp1.*gwtExp;
H2 = jacobExp.*Ttemp2.*gwtExp;

H1 = sum(H1, 4);
H2 = sum(H2, 4);

HVis1 = reshape(H1, D, []);
HVis2 = reshape(H2, D, []);

end

function [HVis1, HVis2] = ViscoKernelB2D(r, dr, drdn, const1, const4, numGauss, jacob, gwt, D)
% ViscoKernelB2D: Calculates the viscoelastic kernel matrices HVis1 and HVis2 
% for the non-singular part of the viscoelastic integral.
%
% Inputs:
%   - r, dr, drdn: Kernel parameters.
%   - const1, const4: Elasticity constants.
%   - numGauss: Number of Gauss points.
%   - jacob: Jacobian matrix.
%   - gwt: Gauss weights.
%   - D: Dimensionality of the problem (2 for 2D).
%
% Outputs:
%   - HVis1, HVis2: Viscoelastic matrices H1 and H2.

% Reshape and broadcast matrices for element-wise operations
r = repmat(r, [D, 1, 1]);
dr = repmat(dr, [D, 1, 1]);
drdn = repmat(drdn, [D, 1, 1]);
dr = permute(dr, [2, 1, 3]);

% drI [1,1,D,numGauss] (3rd dimension:dr1, dr2)
drA = reshape(dr(:,1,:), [1, 1, D, numGauss]);
drExp = reshape(dr, [D, 1, D, numGauss]);
drTExp = permute(drExp, [2, 1, 3, 4]); % Transpose operation for tensor product
RiRj = drExp .* drTExp; % Tensor product ri * rj

% Identity matrix broadcasted over Gauss points
I = eye(D);
IExp = repmat(I, [1, 1, D, numGauss]);
IAExp = reshape(repmat(I, [1, 1, numGauss]), [D, 1, D, numGauss]);

% Compute ikj tensor and its transpose
ikj = IAExp .* drTExp;
ikjT = permute(ikj, [2, 1, 3, 4]);

% Compute coefficient based on dimension
if D == 3
    coeff = drdn ./ r;
elseif D == 2
    coeff = drdn; % For 2D, coefficient doesn't involve r
else
    error('Wrong Dimension! Only D=2 or D=3 are supported.')
end

coeffExp = permute(coeff, [2, 1, 3]); % Adjust dimensions for broadcasting
coeffExp = reshape(coeffExp, [1, 1, D, numGauss]);

% Broadcast Jacobian and Gauss weights for integration
jacobExp = reshape(jacob, [1, 1, 1, numGauss]);
jacobExp = repmat(jacobExp, [1, 1, D, 1]);

gwtExp = reshape(gwt(:), [1, 1, 1, numGauss]);
gwtExp = repmat(gwtExp, [1, 1, D, 1]);

% Compute the viscoelastic kernel matrices HVis1 and HVis2
Ttemp1 = -const1 * coeffExp .* (const4 * (ikj + drA .* IExp) - ikjT + D * drA .* RiRj);
Ttemp2 = -2 * const4 * const1 * coeffExp .* ikjT;

HVis1 = sum(jacobExp .* Ttemp1 .* gwtExp, 4);
HVis2 = sum(jacobExp .* Ttemp2 .* gwtExp, 4);

% Reshape the results into final output matrices
HVis1 = reshape(HVis1, D, []);
HVis2 = reshape(HVis2, D, []);

end

function RadialInte = computeCollinear(RNorm, s, dA, dr, R, D, a, b)
    % Case 2.2.1: RBar + s = 0 (A, p and Q are all in a line).
    ba1 = (b - a)'; 
    ba2 = ba1 .* (b + a)';
    ba3 = (b.^3 - a.^3)';
    ba4 = ba2 .* (b.^2 + a.^2)';
    
    r1Abs = abs(a' - RNorm);
    r2Abs = abs(b' - RNorm);
    
    dr = reshape(dr, D, []);
    R = reshape(R, D, []);
    
    r213 = (r2Abs.^3 - r1Abs.^3);
    baAbs = ((b' - RNorm).*r2Abs - (a' - RNorm).*r1Abs);
    
    RadialInte = 12/dA^4 * (-dr .* ba4 / 4 ...
              - (R + 2*s.*dr).* ba3 / 3 ...
              - (RNorm.^2.*dr + 2*s.*R + dA^2.*dr).* ba2/2 ...
              - (RNorm.^2 + dA^2).* R .* ba1 ...
              + 2*dA*dr.* r213 / 3 ...
              + dA * (R - s.*dr) .* baAbs);
end

function RadialInte = computeNonCollinear(RBarNorm, s, dA, dr, RBar, D, a, b)
    % Case 2.2.2: A, p and Q are not in a line.
    aSqua = a.^2;
    bSqua = b.^2;

    R1 = sqrt(aSqua' + 2 * s .* a' + RBarNorm.^2);
    R2 = sqrt(bSqua' + 2 * s .* b' + RBarNorm.^2);
    dr = reshape(dr, D, []);
    RBar = reshape(RBar, D, []);

    ba1 = (b - a)'; 
    ba2 = (bSqua - aSqua)';  
    ba3 = (b.^3 - a.^3)';  
    ba4 = ((bSqua - aSqua) .* (bSqua + aSqua))';   
    R213 = (R2.^3 - R1.^3);             
    b2a1 = (b' + s).* R2 - (a' + s).* R1; 

    logTerm = (RBarNorm.^2 - s.^2) .* log(abs((R2 + b' + s) ./ (R1 + a' + s))); 

    RadialInte = 12 / dA^4 * (-dr / 4 .* ba4 ...
              - (RBar + 2 * s .* dr) .* ba3 / 3 ...
              - (RBarNorm.^2 .* dr + 2 * s .* RBar + dA^2 * dr) .* ba2 / 2 ...
              - (RBarNorm.^2 + dA^2) .* RBar .* ba1 ...
              + 2 / 3 * dA * dr .* R213 ...
              + dA * (RBar - s .* dr) .* (b2a1 + logTerm));
end