function [H, G] = elasticMatrixInsideRVE2D(nurbsStr, numR, innerPts)
% Assemble System Matrices
%
% Inputs:
%   nurbsStr: Structure of boundary structures containing
%       geometric and material information.
%   numR: Number of regular integration points.
%   inner_points: The points inside.
%
% Outputs:
%   H: Elastic matrix for U.
%   G: Elastic matrix for T.

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
NumIn = size(innerPts,1);
H = zeros(NumIn * D, numU * D);
G = zeros(NumIn * D, numT * D);

% Gaussian Integration
[ptsR, weightR] = lgwt(numR, -1, 1);

for i = 1:NumIn
    row = i*D-D+1 : i*D;
    sourcePt = innerPts(i, :);

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

        % Regular Integration
        [HSub, GSub] = integralRegular2D(elemCoords, gloIdx, sourcePt, ...
            range, p, knotU, weights, const1, const2, const3, const4, ptsR, weightR);

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
    end

    if i == 2
        tEnd = toc(time);
        fprintf(['------------------------------------------------\n' 'Elastic matrix for inner\n' ...
            'One loop time: %f4 s\n' 'Loop number: %d\n' 'Total estimated time %d s = %f2 min\n'], ...
            tEnd, NumIn, NumIn * tEnd, (NumIn * tEnd / 60));
    end
end

end
