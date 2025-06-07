function [HVis1, HVis2] = viscoelasticMatrix3D(nurbsStr, numR, numS, appliedPts, dA, collocNew)
% Function: viscoelasticMatrix3D
% Description: Assembles the viscoelastic system matrices for a 3D problem.
%
% Inputs:
%   nurbsStr: Structure of boundary structures containing geometric and
%             material information.
%   numR: Number of regular integration points.
%   numS: Number of singular integration points.
%   appliedPts: Applied points in RIM, containing boundary points and inner points.
%   dA: Support domain size in RIM.
%   collocNew: new collocation pts.
%
% Outputs:
%   HVis1: Viscoelastic matrix for the deviatoric part (ε*kij × εij).
%   HVis2: Viscoelastic matrix for the volumetric part (ε*kii × εjj).

% Constants
D = 3; % Dimension

% Material Parameters
nu = nurbsStr(1).nu;
E = nurbsStr(1).E;
mu = E / (2 * (1 + nu));
const1 = 1 / (8 * (D - 1) * pi * (1 - nu) * mu);
const2 = 3 - 4 * nu;
const3 = 1 / (4 * (D - 1) * pi * (1 - nu));
const4 = 1 - 2 * nu;

% Initialize indices
[indU, ~, indNumU, indNumT, indNumB] = initializeIndices(nurbsStr, collocNew);
numU = indNumU(end);
numT = indNumT(end);
indElem = [];
numElem = 0;
for i = 1:length(nurbsStr)
    tempE = nurbsStr(i).numElem;
    [indElem] = [indElem; ones(tempE, 1) * i, (1:tempE)'];
    numElem = nurbsStr(i).numElem + numElem;
end

% Initialize Matrices
numApplied = size(appliedPts, 1);
HVis1 = zeros(numApplied * D, (numApplied + D + 1) * D);
HVis2 = zeros(numApplied * D, (numApplied + D + 1) * D);

% Gaussian Integration
[ptsR, weightR] = lgwt(numR, -1, 1);
[ptsS, weightS] = lgwt(numS, -1, 1);
[gptUR, gptVR] = meshgrid(ptsR, ptsR);
[gwtUR, gwtVR] = meshgrid(weightR, weightR);
[gptUS, gptVS] = meshgrid(ptsS, ptsS);
[gwtUS, gwtVS] = meshgrid(weightS, weightS);

for i = 1:numApplied

    row = i * D - D + 1:i * D;
    sourcePt = appliedPts(i, :);
    ip = 0;
    xiParam = 0;

    if i <= numU        % On the boundary.
        ip = indNumB(indU(i, 1)) + indU(i, 2);
        c = indU(i, 3);
        xiParam = nurbsStr(ip).collocXi(c, :);
        bouRange = [indNumB(indU(i, 1))+1, indNumB(indU(i, 1)+1)];
    end

    if i == 2
        time = tic;
    end

    for e = 1:numElem
        id = indElem(e, 1);
        el = indElem(e, 2);
        range = nurbsStr(id).elemRange(el, :);
        gloIdx = nurbsStr(id).globIdx(el, :);

        % Element Coordinates
        elemCoords = nurbsStr(id).ctrlPts(gloIdx, 1:D);
        p = nurbsStr(id).p;
        q = nurbsStr(id).q;
        knotU = nurbsStr(id).knotU;
        knotV = nurbsStr(id).knotV;
        weights = nurbsStr(id).weights;
        uConn = nurbsStr(id).uConn(el, :);

        % Singular Check
        if i <= numU
            flagLine = ismember(nurbsStr(ip).globU(c), nurbsStr(id).globU(uConn))&& id ~= ip;
            if flagLine
                indTemp =  uConn(nurbsStr(id).globU(uConn) == nurbsStr(ip).globU(c));
                xiLine(1) = nurbsStr(id).collocXi(indTemp,1);
                xiLine(2) = nurbsStr(id).collocXi(indTemp,2);
            else
                xiLine = [0,0];
            end

            % Judge the type of integral
            [xi_u, xi_v, flag] = JudgeSingular3D(xiParam, range, ip, id, bouRange, flagLine, xiLine);
        else
            flag = 0;
            xi_u = 0;
            xi_v = 0;
        end

        [HVis1Sub1, HVis2Sub1, HVis1Sub2, HVis2Sub2] = integralViscoelastic3D(elemCoords, ...
            gloIdx, sourcePt, xi_u, xi_v, range, p, q, knotU, knotV, weights, const1, const2, const3, ...
            const4, gptUR, gptVR, gptUS, gptVS, gwtUR, gwtVR, gwtUS, gwtVS, appliedPts, dA, flag, D);

        % Accumulate Matrices
        HVis1(row, :) = HVis1(row, :) + [HVis1Sub1 zeros(D) HVis1Sub2];
        HVis2(row, :) = HVis2(row, :) + [HVis2Sub1 zeros(D) HVis2Sub2];
    end

    if i == 2
        tEnd = toc(time);
        fprintf(['------------------------------------------------\n' 'Viscoelastic matrix\n' ...
            'One loop time: %f4 s\n' 'Loop number: %d\n' 'Total estimated time %d s = %f2 min\n'], ...
            tEnd, numApplied, numApplied*tEnd, (numApplied*tEnd/60));
    end
end

end
