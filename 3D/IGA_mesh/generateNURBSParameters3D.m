function nurbsParams = generateNURBSParameters3D(shapeType, shapeParams)
% Generates NURBS parameters for 3D shapes based on boundary information.
%
% Input:
%   shapeType : Type of shape.
%   shapeParams : Information defining the boundary shape.
%
% Output:
%   nurbsParams (struct): Struct containing NURBS parameters.
%           ______
%          /     /| →3
%         /  5  / |
%        /_____/ 2|
%     4← |     |  /
%        |  1  | /
%        |_____|/
%           ↑
%           6

nurbsParams = struct;

switch shapeType
    case 'cub'
        numPatch = 6;
        lengthCube = shapeParams(4 : 6);
        centerCube = shapeParams(1 : 3);
        [knotU, knotV] = deal(cell(numPatch, 1));
        for i=1:numPatch
            [knotU{i}, knotV{i}] = deal([0, 0, 0, 1, 1, 1]);
        end
        [controlPoints, weights] = deal(cell(numPatch, 1));
        patchLink = [6 2 5 4; 6 3 5 1; 6 4 5 2; 6 1 5 3; 1 2 3 4; 3 2 1 4];
        for j = 1:numPatch
            controlPoints{j} = getPointsCube(lengthCube, centerCube, j);
            weights{j} = ones(9, 1);
            nurbsParams(j).patchLink = patchLink(j, :);
        end
        isClosedU = 0;
        isClosedV = 0;
        p = 2;
        q = 2;
    case 'sph6'
        numPatch = 6;
        r = shapeParams(1);
        o = shapeParams(2 : 4);
        [controlPoints, weights, knotUTemp, knotVTemp] = getPointsSphere6(r, o);
        [knotU, knotV] = deal(cell(numPatch, 1));
        for i=1:numPatch
            [knotU{i}, knotV{i}] = deal(knotUTemp, knotVTemp);
        end
        patchLink = [6 2 5 4; 6 3 5 1; 6 4 5 2; 6 1 5 3; 1 2 3 4; 3 2 1 4];
        for j = 1:numPatch
            nurbsParams(j).patchLink = patchLink(j, :);
        end
        isClosedU = 0;
        isClosedV = 0;
        p = 4;
        q = 4;
    case 'sphP'
        numPatch = 1;
        r = shapeParams(1);
        o = shapeParams(2 : 4);
        knotU = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4];
        knotV = [0, 0, 0, 1, 1, 2, 2, 2];
        controlPoints{1} = getPointsSphereP(r, o);
        sq2 = sqrt(2) / 2;
        weights{1} = [1; sq2; 1; sq2; 1; sq2; 1; sq2; 1; sq2; 0.5; sq2; 0.5; sq2; 0.5; sq2; 0.5; sq2; ...
            1; sq2; 1; sq2; 1; sq2; 1; sq2; 1; sq2; 0.5; sq2; 0.5; sq2; 0.5; sq2; 0.5; sq2; 1; sq2; 1; sq2; 1; sq2; 1; sq2; 1];
        patchLink = [];
        nurbsParams.patchLink = patchLink;
        isClosedU = 1;
        isClosedV = 0;
        p = 2;
        q = 2;
    case 'star_e3'
        numPatch = 6;
        paraStar = shapeParams(4 : 8);
        o = shapeParams(1 : 3);
        [controlPoints, weights, knotU, knotV] = getPointsStarExa3(paraStar, o);
        patchLink = [6 2 5 4; 6 3 5 1; 6 4 5 2; 6 1 5 3; 1 2 3 4; 3 2 1 4];
        for j = 1:numPatch
            nurbsParams(j).patchLink = patchLink(j, :);
        end
        isClosedU = 0;
        isClosedV = 0;
        p = 2;
        q = 2;
    case 'star_e4'
        numPatch = 6;
        paraStar = shapeParams(4 : 8);
        o = shapeParams(1 : 3);
        [controlPoints, weights, knotU, knotV] = getPointsStarExa4(paraStar, o);
        patchLink = [6 2 5 4; 6 3 5 1; 6 4 5 2; 6 1 5 3; 1 2 3 4; 3 2 1 4];
        for j = 1:numPatch
            nurbsParams(j).patchLink = patchLink(j, :);
        end
        isClosedU = 0;
        isClosedV = 0;
        p = 2;
        q = 2;
    case 'cyli_e4'
        numPatch = 6;
        paraStar = shapeParams(4 : 6);
        o = shapeParams(1 : 3);
        [controlPoints, weights, knotU, knotV] = getPointsCylinderExa4(paraStar, o);
        patchLink = [6 2 5 4; 6 3 5 1; 6 4 5 2; 6 1 5 3; 1 2 3 4; 3 2 1 4];
        for j = 1:numPatch
            nurbsParams(j).patchLink = patchLink(j, :);
        end
        isClosedU = 0;
        isClosedV = 0;
        p = 2;
        q = 2;
    otherwise
        error('No such case!')
end

% Create struct to store NURBS parameters
for j = 1:numPatch
    nurbsParams(j).ctrlPts = controlPoints{j};
    nurbsParams(j).knotU = knotU{j};
    nurbsParams(j).knotV = knotV{j};
    nurbsParams(j).weights = weights{j};
    nurbsParams(j).isClosedU = isClosedU;
    nurbsParams(j).isClosedV = isClosedV;
    nurbsParams(j).p = p;
    nurbsParams(j).q = q;
end

end
