function avgStress = getAverageStress(nurbsStr, trac)
% Constants
D = 2; % Dimension

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

trac = [trac(1:2:end), trac(2:2:end)];

% Initialize Matrices
avgStress = zeros(1, 3);

% Gaussian Integration
numG = 12;
[ptsG, weightG] = lgwt(numG, -1, 1);

% Out boundary.
for e = 1:nurbsStr(1).numElem
    id = 1;
    el = e;
    range = nurbsStr(id).elemRange(el, :);
    gloIdx = nurbsStr(id).globIdx(el, :);
    elemCoords = nurbsStr(id).ctrlPts(gloIdx, 1:D);
    uConn = nurbsStr(id).uConn(el, :);
    tConn = nurbsStr(id).tConn(el, :);

    % Element Coordinates
    ctrlPts = nurbsStr(id).ctrlPts;
    p = nurbsStr(id).p;
    knotU = nurbsStr(id).knotU;
    weights = nurbsStr(id).weights;

    % Ti and Xj
    Ti = trac(tConn, 1:D);
    Xj = ctrlPts(uConn, 1:D);

    % Integration
    avgStressSub = integralAvgStress2D(elemCoords, gloIdx, range, p, ...
        knotU, weights, ptsG, weightG, Ti, Xj);

    % Accumulate Matrices
    avgStress = avgStress + avgStressSub;
end


end