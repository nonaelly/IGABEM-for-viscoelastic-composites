function [H, G] = elasticMatrixBoundary3D(nurbsStr, collocNew, ...
    numR, numS, tranU, appliedPts)
% Assemble System Matrices
%   A viscoelastic matrix with several inclusions.
%
% Inputs:
%   nurbsStr: Structure of boundary structures containing
%       geometric and material information.
%   collocNew: Struct containing new collocation points and global indices.
%   numR: Number of regular integration points.
%   numS: Number of singular integration points.
%   tranU: Transformation matrix for displacement.
%   appliedPts: Applied points in radial integration method.
%
% Outputs:
%   H: Elastic matrix for U.
%   G: Elastic matrix for T.

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

% Initialize
[indU, ~, indNumU, indNumT, indNumB] = initializeIndices(nurbsStr, collocNew);
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

% Gaussian Integration
[ptsR, weightR] = lgwt(numR, -1, 1);
[gptUR, gptVR] = meshgrid(ptsR, ptsR);
[gwtUR, gwtVR] = meshgrid(weightR, weightR);

for i = 1:numU
    row = i * D - D + 1 : i * D;
    ip = indNumB(indU(i, 1)) + indU(i, 2);
    c = indU(i, 3);
    xiParam = nurbsStr(ip).collocXi(c, :);
    bouRange = [indNumB(indU(i, 1))+1, indNumB(indU(i, 1)+1)];
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
        q = nurbsStr(id).q;
        knotU = nurbsStr(id).knotU;
        knotV = nurbsStr(id).knotV;
        weights = nurbsStr(id).weights;

        flagLine = ismember(nurbsStr(ip).globU(c), nurbsStr(id).globU(uConn)) && id ~= ip;
        if flagLine
            indTemp =  uConn(nurbsStr(id).globU(uConn) == nurbsStr(ip).globU(c));
            xiLine(1) = nurbsStr(id).collocXi(indTemp,1);
            xiLine(2) = nurbsStr(id).collocXi(indTemp,2);
        else
            xiLine = [0,0];
        end
        
        % Judge the type of integral
        [xi_u, xi_v, flag] = JudgeSingular3D(xiParam, range, ip, id, bouRange, flagLine, xiLine);

        if flag
            % Singular Integration (Power series expression).
            [xi, eta] = getXiParamter(xi_u, xi_v, range);
            
            [HSubO, GSubO] = integralSingular3D(numS, xi, eta, p, q, range, ...
                elemCoords, collocPt, knotU, knotV, weights, const1, const2, const3, const4, id);
        else
            % Regular Integration
            [HSubO, GSubO] = integralRegular3D(elemCoords, collocPt, range, p, q, knotU, knotV, weights, ...
                const1, const2, const3, const4, gptUR, gptVR, gwtUR, gwtVR, D, id);
        end

        uConnBou = nurbsStr(id).globU(uConn) + indNumU(nurbsStr(id).bouInd);
        tConnBou = nurbsStr(id).globT(tConn) + indNumT(nurbsStr(id).bouInd);
        colU = zeros(1,D*(p+1)*(q+1));
        colT = zeros(1,D*(p+1)*(q+1));
        for k=1:D
            colU(k:D:end) = uConnBou*D-D+k;
            colT(k:D:end) = tConnBou*D-D+k;
        end

        % Accumulate Matrices
        H(row, colU) = H(row, colU) + HSubO;
        G(row, colT) = G(row, colT) + GSubO;
    end

    if i == 2
        tEnd = toc(time);
        fprintf(['------------------------------------------------\n' 'Elastic matrix for boundary\n' ...
            'One loop time: %f4 s\n' 'Loop number: %d\n' 'Total estimated time %d s = %f2 min\n'], ...
            tEnd, numU, numU * tEnd, (numU * tEnd / 60));
    end
end

% Rigid body method for IGABEM
H = H / (tranU);

for i = 1:numU
    H((i * 3 - 2):(i * 3), (i * 3 - 2):(i * 3)) = 0; % Set the diagonal (3x3) to 0
end
for i = 1:numU
    rowI = i * 3 - 2:i * 3;
    Hbii = -[sum(H(rowI, 1:3:end), 2), sum(H(rowI, 2:3:end), 2), ...
        sum(H(rowI, 3:3:end), 2)];
    H(rowI, rowI) = Hbii;
end
H = H * (tranU);

end
