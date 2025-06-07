function meshParams = generateMesh3D(nurbs, lamda)
% Generates mesh parameters for a 3D NURBS shape.
%
% Input:
%   nurbs (struct): NURBS parameters.
%   shapeType (string): Type of the shape ('cube', 'sphereP', etc.).
%   lamda (scalar): Discontinuous element parameter.
%
% Output:
%   meshParams (struct): Mesh parameters including control points, weights, and connectivity.

% Extract NURBS parameters
knotU = nurbs.knotU;
knotV = nurbs.knotV;
ctrlPts = nurbs.ctrlPts;
weights = nurbs.weights;
p = nurbs.p;
q = nurbs.q;
refinement = nurbs.refinement;
isClosedU = nurbs.isClosedU;
isClosedV = nurbs.isClosedV;
patchLink = nurbs.patchLink;

% Normalize knot vectors
knotU = knotU / max(knotU);
knotV = knotV / max(knotV);

% Initialize mesh parameters
meshParams = struct;

% h-Refinement
[ctrlPts, knotU, knotV, weights] = hRefinement3D(p, q, refinement, knotU, knotV, ctrlPts, weights);

% Get connectivity matrices
[globIdx, uConn, tConn, t2uConn, elemRange] = connectMatrix3D(p, q, knotU, knotV, isClosedU, isClosedV);

% Collocation points
m = length(knotU) - p - 1;
n = length(knotV) - q - 1;
collocXi = zeros((m-isClosedU)*(n-isClosedV), 2);
collocPts = zeros((m-isClosedU)*(n-isClosedV), 3);

% Compute collocation points
tole = 1e-8;
for i = 1 : n-isClosedV
    for j = 1 : m-isClosedU
        ind = j + (i-1)*(m-isClosedU);
        xi = sum(knotU((j + 1):(j + p))) / p;
        eta = sum(knotV((i + 1):(i + q))) / q;
        collocXi(ind,:) = [xi, eta];
        % Discontinuous element.
        for k=1:2
            if (abs(collocXi(ind,k))<tole || abs(1-collocXi(ind,k))<tole)
                collocXi(ind,k) = abs(1-lamda-collocXi(ind,k));
            end
        end
        collocPts(ind, :) = NURBSinterpolation3d(collocXi(ind,1), collocXi(ind,2), p, q, knotU, knotV, ctrlPts, weights');
    end
end

numCollocU = size(collocPts,1);
numCollocT = length(t2uConn);
numElemU = length(unique(knotU)) - 1;
numElemV = length(unique(knotV)) - 1;
numElem = numElemU * numElemV;

normCollocT = zeros(size(collocPts, 1), 3);
for i = 1 : size(collocPts, 1)
    xi = collocXi(i,:);
    for e = 1 : numElem
        range = elemRange(e, :);
        gloIdx = globIdx(e, :);
        elemCoords = ctrlPts(gloIdx, 1:3);
        if range(2) >= xi(1) && xi(1) >=range(1) && range(4) >= xi(2) && xi(2) >= range(3)
            [Xi_u, Xi_v] = getXiParamter(xi(1), xi(2), range);
            [~, dNu, dNv] = NURBS2DBasisDers([Xi_u, Xi_v], p, q, knotU, knotV, weights');
            dxydxi1 = dNu * elemCoords;
            dxydxi2 = dNv * elemCoords;
            normals = cross(dxydxi1, dxydxi2);
            normals = normals / norm(normals);
            normCollocT(i, :) = normals;
            break
        end
    end
end

% Structure
meshParams.ctrlPts = ctrlPts;
meshParams.weights = weights;
meshParams.knotU = knotU;
meshParams.knotV = knotV;
meshParams.collocPts = collocPts;
meshParams.collocXi = collocXi;
meshParams.p = p;
meshParams.q = q;
meshParams.refinement = refinement;
meshParams.elemRange = elemRange;
meshParams.globIdx = globIdx;
meshParams.uConn = uConn;
meshParams.tConn = tConn;
meshParams.t2uConn = t2uConn;
meshParams.patchLink = patchLink;
meshParams.numCollocU = numCollocU;
meshParams.numCollocT = numCollocT;
meshParams.numElemU = numElemU;
meshParams.numElemV = numElemV;
meshParams.numElem = numElem;
meshParams.normCollocT = normCollocT;

end
