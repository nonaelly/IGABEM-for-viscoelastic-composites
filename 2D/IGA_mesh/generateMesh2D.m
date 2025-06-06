function meshParameters = generateMesh2D(nurbsParameters)
% Generates the mesh parameters for a 2D problem using NURBS.
%
% Input:
%   nurbsParameters (struct): Struct containing parameters for NURBS curves.
%
% Output:
%   meshParameters (struct): Struct containing mesh parameters for the 2D problem.

% Initialize variables
knotVector = nurbsParameters.knotVector;
knotVector = knotVector / max(knotVector);
controlPoints = nurbsParameters.controlPoints;
weights = nurbsParameters.weights;
p = nurbsParameters.orderNURBS;
refinement = nurbsParameters.refinement;
isClosed = nurbsParameters.isClosed;

% h-Refinement
[controlPoints, knotVector, weights] = hRefinement(p, refinement, ...
    knotVector, controlPoints, weights);

% Get connection matrix
[globalBasisIndex, UConnection, TConnection, elementRange] = connectMatrix(p, ...
    knotVector, controlPoints, weights, isClosed);

% Get the connection between U and T
numElement = size(globalBasisIndex, 1); % Number of elements
T2UConnection = zeros(max(TConnection(end, :)), 1);
for e = 1:numElement
    T2UConnection(TConnection(e, :)) = UConnection(e, :);
end

% Get the collocation points in parameter and physical space
numControlPoints = size(controlPoints, 1);
collocationXi = zeros(numControlPoints - isClosed, 1);
collocationPoints = zeros(numControlPoints - isClosed, 2);

% Parameter
for i = 1:numControlPoints - isClosed
    collocationXi(i) = sum(knotVector((i + 1):(i + p))) / p;
end

% Physical
for i = 1:(numControlPoints - isClosed)
    collocationPoints(i, 1) = NURBSinterpolation(collocationXi(i), p, knotVector, ...
        controlPoints(:, 1)', weights);
    collocationPoints(i, 2) = NURBSinterpolation(collocationXi(i), p, knotVector, ...
        controlPoints(:, 2)', weights);
end

% Basic information for U and T on collocation points
numCollocationU = size(collocationPoints, 1);
numCollocationT = length(T2UConnection);

% Norm vectors of every collocation point (2 vectors)
normsOnCollocation = getNormCollocationPoints(elementRange, globalBasisIndex, ...
    collocationXi, UConnection, controlPoints, weights, knotVector, p, isClosed);

% Norm vectors of every traction point (1 vector)
normCollocationT = zeros(numCollocationT, 2);
e = 1;
for i = 1:numCollocationT
    j = T2UConnection(i);
    if i == TConnection(e, end)
        normCollocationT(i, :) = normsOnCollocation(j, 1:2);
        e = e + 1;
    else
        normCollocationT(i, :) = normsOnCollocation(j, 3:4);
    end
end

% Create struct to store mesh parameters with shorter variable names
meshParameters = struct;
meshParameters.ctrlPts = controlPoints;
meshParameters.weights = weights;
meshParameters.knotU = knotVector;
meshParameters.collocPts = collocationPoints;
meshParameters.collocXi = collocationXi;
meshParameters.p = p;
meshParameters.refine = refinement;
meshParameters.elemRange = elementRange;
meshParameters.globIdx = globalBasisIndex;
meshParameters.uConn = UConnection;
meshParameters.tConn = TConnection;
meshParameters.t2uConn = T2UConnection;
meshParameters.numCollocU = numCollocationU;
meshParameters.numCollocT = numCollocationT;
meshParameters.numElem = numElement;
meshParameters.normsColloc = normsOnCollocation;
meshParameters.normCollocT = normCollocationT;
meshParameters.isClosed = isClosed;
end
