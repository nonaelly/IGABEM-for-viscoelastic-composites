function [ctrlPts, weights, knotU, knotV] = getPointsStarExa3(paraStar, o)
[H1, H2, H3, R1, R2] = deal(paraStar(1), paraStar(2), paraStar(3), paraStar(4), paraStar(5));

addpath ../Pre
iS2 = 1/sqrt(2);
RMid = R1 - H2;
knotsForInsertion = [1/5 1/5 2/5 2/5 3/5 3/5 4/5 4/5] ;

% Outer circle
nodeSeg = [0,R1;R1,0;];
nurbsTemp = generateNURBSParameters(1, 1, nodeSeg, o, 2, 0);
ctrlPtsTemp = nurbsTemp.controlPoints;
weightsOut = nurbsTemp.weights;
ctrlPtsTemp = [ctrlPtsTemp, zeros(size(ctrlPtsTemp, 1), 1)];
knotTemp = [0,0,0,1,1,1];
pTemp = 2;
nurbsOut = NURBS1D(knotTemp, pTemp, ctrlPtsTemp, weightsOut);

% Outer circle refine
nurbsOut.p_refine_curve(0);
nurbsOut.h_refine_curve(knotsForInsertion);

numPatch = 6;
weights = cell(numPatch, 1);
weights{1} = [1; iS2; 1; 1; 1; iS2; 1; 1; 1; iS2; 1; 1; iS2; 1; 1; 1; iS2; 1; 1; 1; iS2; 1;];
weights{3} = [nurbsOut.weights; nurbsOut.weights;];
weights{2} = ones(4,1);
weights{4} = ones(4,1);
weights{5} = [1; iS2; 1; 1; 1; iS2; 1; 1; 1; iS2; 1; nurbsOut.weights];
weights{6} = [nurbsOut.weights; 1; iS2; 1; 1; 1; iS2; 1; 1; 1; iS2; 1];

ctrlPts = cell(numPatch,1);

ctrlPtsTemp = [0,RMid; R2,RMid; R2,RMid-R2; R2,(RMid+R2)/2; R2,2*R2; R2,R2; 2*R2,R2; ...
    (RMid+R2)/2,R2; RMid-R2,R2; RMid,R2; RMid,0];
ctrlPts{1}=[ctrlPtsTemp, zeros(size(ctrlPtsTemp, 1), 1); ctrlPtsTemp, H3*ones(size(ctrlPtsTemp, 1), 1)];

nodeSeg = [RMid,0; R1,0];
nurbsTemp = generateNURBSParameters(1, 0, nodeSeg, [], 1, 0);
ctrlPtsTemp = nurbsTemp.controlPoints;
ctrlPts{2}=[ctrlPtsTemp, zeros(size(ctrlPtsTemp, 1), 1); ctrlPtsTemp, H3*ones(size(ctrlPtsTemp, 1), 1)];

ctrlPtsTemp1 = nurbsOut.controlPts;
ctrlPtsTemp1(:,3) = ctrlPtsTemp1(:,3) + H3;
ctrlPts{3}=[flip(nurbsOut.controlPts); flip(ctrlPtsTemp1); ];

nodeSeg = [0,R1; 0,RMid];
nurbsTemp = generateNURBSParameters(1, 0, nodeSeg, [], 1, 0);
ctrlPtsTemp = nurbsTemp.controlPoints;
ctrlPts{4}=[ctrlPtsTemp, zeros(size(ctrlPtsTemp, 1), 1); ctrlPtsTemp, H3*ones(size(ctrlPtsTemp, 1), 1)];

ctrlPtsTemp = [0,RMid; R2,RMid; R2,RMid-R2; R2,(RMid+R2)/2; R2,2*R2; R2,R2; 2*R2,R2; ...
    (RMid+R2)/2,R2; RMid-R2,R2; RMid,R2; RMid,0];
ctrlPts{5}=[ctrlPtsTemp(:,1) ctrlPtsTemp(:,2) H3*ones(size(ctrlPtsTemp,1),1); ctrlPtsTemp1];

ctrlPts{6}=[nurbsOut.controlPts; ctrlPtsTemp(:,1) ctrlPtsTemp(:,2) zeros(size(ctrlPtsTemp,1),1)];

% genenrate each NURBS patch information object
q = 1;
vKnot = [0 0 1 1];
shellNurbs = NURBS(numPatch);
for i = 1:numPatch
    if ismember(i, [2,4])
        uKnot = [0 0 1 1];
        p = 1;
    else
        uKnot = [0 0 0 knotsForInsertion 1 1 1];
        p = 2;
    end
    shellNurbs(i) = NURBS(uKnot, vKnot, p, q, ctrlPts{i}, weights{i});
end

% p-refinement
pAdd = 0;
pAdd1 = pAdd + 1;
qAdd = 1;
for i = 1 : numPatch
    if ismember(i, [2,4])
        shellNurbs(i).p_refine_surface(pAdd1,qAdd);
    else
        shellNurbs(i).p_refine_surface(pAdd,qAdd);
    end
end

% h-refinement
nx = 25;
knotx = [];
for i = 1 : 4
    [knotx] = [knotx, (i:5:nx-1)/nx];
end
knotx = sort(knotx);
ny = 6;
nz = 4;
knoty = (1:ny-1)/ny;
knotz = (1:nz-1)/nz;

knotx = [];
knoty = [];
knotz = [];

insertKnot = {knotx,knotz; knoty,knotz; knotx,knotz; knoty,knotz; knotx,knoty; knotx,knoty};
for i = 1:numPatch
    shellNurbs(i).h_refine_surface(insertKnot{i,1}, insertKnot{i,2});
end

ctrlPts = cell(numPatch,1);
weights = cell(numPatch, 1);
[knotU, knotV] = deal(cell(numPatch, 1));
for i = 1:numPatch
    knotU{i} = shellNurbs(i).uKnot;
    knotV{i} = shellNurbs(i).vKnot;
    ctrlPts{i} = shellNurbs(i).controlPts;
    weights{i} = shellNurbs(i).weights;
end
rmpath ../Pre

end