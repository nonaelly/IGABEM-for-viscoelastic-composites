function nurbsStr = getStructure3D(geomData, E, nu, refinement, lamda, numPlot)
% Generates the struct for a 3D problem and plots the model photo.
%
% Input:
%   geomData (cell array): Geometry information.
%   p (vector): NURBS order for each boundary segment.
%   E (vector): Young's modulus for each region.
%   nu (vector): Poisson's ratio for each region.
%   refinement (vector): Order of refinement for each boundary segment.
%   lamda (scalar): Discontinuous element parameter.
%   numPlot (scalar): Figure number for plotting.
%
% Output:
%   nurbsStr (struct array): Array of structs containing information for each boundary.
%       Each struct includes data for BEM mesh generation and plotting.

% Initialize variables
numBoundary = length(geomData);
patchInd = 1;

% Loop over each boundary
for i = 1:numBoundary
    shapeType = geomData{i}{1};
    shapeParams = geomData{i}{2};
    nurbs = generateNURBSParameters3D(shapeType,shapeParams);

    % Generate mesh
    for j=1:length(nurbs)
        switch j
            case {1, 3}
                nurbs(j).refinement = refinement(i, [1, 3]);
            case {2, 4}
                nurbs(j).refinement = refinement(i, [2, 3]);
            case {5, 6}
                nurbs(j).refinement = refinement(i, [1, 2]);
        end
        boundaryStr(j) = generateMesh3D(nurbs(j), lamda(i));
    end

    % Assign material properties
    for j=1:length(nurbs)
        boundaryStr(j).E = E(i);
        boundaryStr(j).nu = nu(i);
        boundaryStr(j).bouInd  = i;
    end

    % Store in structure array
    for j=1:length(nurbs)
        nurbsStr(patchInd) = boundaryStr(j);
        patchInd = patchInd + 1;
    end

    % Plot the mesh
%     if i>1
    plotMesh3D(boundaryStr,numPlot);
%     end

    clear boundaryStr nurbs
end



end
