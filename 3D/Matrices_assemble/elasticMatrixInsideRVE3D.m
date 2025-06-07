function [H, G] = elasticMatrixInsideRVE3D(nurbsStr, numGauR, innerPoints, collocNew)
% Assemble System Matrices
%   A viscoelastic matrix with several inclusions.
%
% Inputs:
%   nurbsStr: Structure of boundary structures containing
%       geometric and material information.
%   numR: Number of regular integration points.
%   innerPoints: Inner points in radial integration method.
%
% Outputs:
%   H: Elastic matrix for U for the RVE's matrix.
%   G: Elastic matrix for T for the RVE's matrix.


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
[indU, indT, indNumU, indNumT, indNumB] = initializeIndices(nurbsStr, collocNew);
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
numIn = size(innerPoints,1);
H = zeros(numIn * D, numU * D);
G = zeros(numIn * D, numT * D);

% Gaussian Integration
[ptsR, weightR] = lgwt(numGauR, -1, 1);
[gptUR, gptVR] = meshgrid(ptsR, ptsR);
[gwtUR, gwtVR] = meshgrid(weightR, weightR);

for i = 1:numIn
    row = i * D - D + 1:i * D;
    sourcePt = innerPoints(i, :);

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

        % Regular Integration
        [HSub, GSub] = integralRegular3D(elemCoords, sourcePt, range, p, q, knotU, knotV, weights, ...
            const1, const2, const3, const4, gptUR, gptVR, gwtUR, gwtVR, D, id);

        uConnBou = nurbsStr(id).globU(uConn) + indNumU(nurbsStr(id).bouInd);
        tConnBou = nurbsStr(id).globT(tConn) + indNumT(nurbsStr(id).bouInd);
        colU = zeros(1,D*(p+1)*(q+1));
        colT = zeros(1,D*(p+1)*(q+1));
        for k=1:D
            colU(k:D:end) = uConnBou*D-D+k;
            colT(k:D:end) = tConnBou*D-D+k;
        end

        % Accumulate Matrices
        coefficient = 1 * (id < 7) + (1 - nurbsStr(id).E/E) * (id > 6);
        H(row, colU) = H(row, colU) + HSub * coefficient;
        G(row, colT) = G(row, colT) + GSub * (id < 7);
    end

    if i == 2
        tEnd = toc(time);
        fprintf(['------------------------------------------------\n' 'Elastic matrix for inner\n' ...
            'One loop time: %f4 s\n' 'Loop number: %d\n' 'Total estimated time %d s = %f2 min\n'], ...
            tEnd, numIn, numIn * tEnd, (numIn * tEnd / 60));
    end
end

end
