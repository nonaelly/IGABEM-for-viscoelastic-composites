function [indU, indT, indNumU, indNumT, indNumB] = initializeIndices(nurbsStr,collocNew)
numBou  = length(collocNew);
indNumU = zeros(numBou+1,1);
indNumT = zeros(numBou+1,1);
indNumB = zeros(numBou+1,1);
indU = []; % boundary; patch; number.
indT = []; % boundary; patch; number.

% A discontinuous cube with continuous/discontinuous spheres.
for i=1:length(collocNew)

    indNumB(i+1) = indNumB(i) + sum([nurbsStr.bouInd] == i);
    indNumU(i+1) = indNumU(i) + size(collocNew(i).collocPts,1);
    indU = [indU; [i * ones(size(collocNew(i).collocPts,1),1), collocNew(i).infoPts]];

    if i == 1 % Cube, or use an index to identify the geometry.
        tempT = 0;
        infoPts = [];
        cubeStr = nurbsStr([nurbsStr.bouInd] == i);
        for j = 1:length(cubeStr)
            tempInfo = [j* ones(cubeStr(j).numCollocT,1),(1:cubeStr(j).numCollocT)'];
            infoPts = [infoPts; tempInfo];
            tempT = tempT + cubeStr(j).numCollocT;
        end
        indNumT(i+1) = indNumT(i) + tempT;
        indT = [indT; [i * ones(tempT,1), infoPts]];
    else
        indNumT(i+1) = indNumT(i) + size(collocNew(i).collocPts,1);
        indT = [indT; [i * ones(size(collocNew(i).collocPts,1),1), collocNew(i).infoPts]];
    end
end
end