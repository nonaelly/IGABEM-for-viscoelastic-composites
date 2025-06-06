function [H, G, HOri, GOri] = elasticMatrixBoundaryRVE2D(nurbsStr, numR, numS, appliedPts)
% Assemble System Matrices for RVE.
%
% Inputs:
%   nurbsStr: Structure of boundary structures containing
%       geometric and material information.
%   numR: Number of regular integration points.
%   numS: Number of singular integration points.
%   appliedPts: Applied points in radial integration method.
%
% Outputs:
%   H: Elastic matrix for U for the matrix (strong singular).
%   G: Elastic matrix for T for the matrix (weakly singular).

% Constants
D = 2; % Dimension

% Material Parameters
nu = nurbsStr(1).nu;
E = nurbsStr(1).E;
mu = E / (2 * (1 + nu));
const1 = 1 / (8 * (D - 1) * pi * (1 - nu) * mu);
const2 = 3 - 4 * nu;
const3 = 1 / (4 * (D - 1) * pi * (1 - nu));
const4 = 1 - 2 * nu;

% Initialize
[indU, indT, indNumU, indNumT, indNumB] = initializeIndices2D(nurbsStr);
numU = indNumU(end);
numT = indNumT(end);
indElem = [];
numElem = 0;
for i = 1:length(nurbsStr)
    tempE = nurbsStr(i).numElem;
    indElem = [indElem; ones(tempE, 1) * i, (1:tempE)'];
    numElem = nurbsStr(i).numElem + numElem;
end

% Initialize Matrices
H = zeros(numU * D);
G = zeros(numU * D, numT * D);
HOri = zeros(numU * D);
GOri = zeros(numU * D, numT * D);

% Gaussian Integration
[ptsR, weightR] = lgwt(numR, -1, 1);
[ptsS, weightS] = lgwt(numS, -1, 1);

for i = 1:numU
    row = i*D-D+1 : i*D;
    ip = indU(i, 1);
    c = indU(i, 2);
    xiParam = nurbsStr(ip).collocXi(c, :);
    collocPt = appliedPts(i, :);

    if i == 2
        time = tic;
    end

    for e = 1:numElem
        id = indElem(e, 1);
        el = indElem(e, 2);
        range = nurbsStr(id).elemRange(el, :);
        gloIdx = nurbsStr(id).globIdx(el, :);
        uConn = nurbsStr(id).uConn(el, :);
        tConn = nurbsStr(id).tConn(el, :);

        % Element Coordinates
        elemCoords = nurbsStr(id).ctrlPts(gloIdx, 1:D);
        p = nurbsStr(id).p;
        knotU = nurbsStr(id).knotU;
        weights = nurbsStr(id).weights;
        
        % Judge the type of integral
        [xi_u, flag] = JudgeSingular2D(nurbsStr, xiParam, range, ip, id, el);

        if flag
            % Singular Integration.
            n1 = nurbsStr(ip).normsColloc(c,1:2);
            n2 = nurbsStr(ip).normsColloc(c,3:4);
            jumpTerm = calculateJumpTerm(n1, n2, nu);

            % (Subtraction of singularity method)
            [HSubOri] = integralSST2D(elemCoords, gloIdx, collocPt, xi_u, ...
                range, p, knotU, weights, jumpTerm, const1, const2, const3, const4, ptsS, weightS);

            if id > 1
                jumpTerm = jumpTerm * (E + nurbsStr(id).E) / (E - nurbsStr(id).E);
            end
            
            [HSub] = integralSST2D(elemCoords, gloIdx, collocPt, xi_u, ...
                range, p, knotU, weights, jumpTerm, const1, const2, const3, const4, ptsS, weightS);

            [GSub] = integralTelles2D(elemCoords, gloIdx, collocPt, xi_u, ...
                range, p, knotU, weights, const1, const2, const3, const4, ptsS, weightS);

            GSubOri = GSub;

        else
            % Regular Integration
            [HSub, GSub] = integralRegular2D(elemCoords, gloIdx, collocPt, ...
                range, p, knotU, weights, const1, const2, const3, const4, ptsR, weightR);
            HSubOri = HSub;
            GSubOri = GSub;
        end

        if id > 1 % Inclusions
            HSub = HSub * (1 - nurbsStr(id).E/E);
            GSub = GSub * 0;
        end

        uConnBou = uConn + indNumU(id);
        tConnBou = tConn + indNumT(id);
        colU = zeros(1,D*(p+1));
        colT = zeros(1,D*(p+1));
        for k=1:D
            colU(k:D:end) = uConnBou*D-D+k;
            colT(k:D:end) = tConnBou*D-D+k;
        end

        % Accumulate Matrices
        H(row, colU) = H(row, colU) + HSub;
        G(row, colT) = G(row, colT) + GSub;
        HOri(row, colU) = HOri(row, colU) + HSubOri;
        GOri(row, colT) = GOri(row, colT) + GSubOri;
    end

    if i == 2
        tEnd = toc(time);
        fprintf(['------------------------------------------------\n' 'Elastic matrix for boundary\n' ...
            'One loop time: %f4 s\n' 'Loop number: %d\n' 'Total estimated time %d s = %f2 min\n'], ...
            tEnd, numU, numU * tEnd, (numU * tEnd / 60));
    end
end

end
