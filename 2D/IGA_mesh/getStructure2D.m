function nurbsStr = getStructure2D(geomData, p, E, nu, refinement, isClosed, numPlot)
% Generates the struct for a 2D problem and plots the model photo.
%
% Input:
%   geomData (struct array): Geometry information.
%   p (vector): NURBS order for each boundary segment.
%   E (scalar): Young's modulus.
%   nu (scalar): Poisson's ratio.
%   refinement (vector): Order of refinement for each boundary segment.
%   isClosed (vector): Indicates whether the curve is closed or open.
%   numPlot (scalar): Number of control points.
%
% Output:
%   nurbsStr (struct array): Array of structs containing information for each boundary.
%       Each struct includes data for BEM mesh generation and plotting.

% Initialize variables
serialSegment = 0; 
numBoundary = length(geomData);

% Loop over each boundary
for i = 1:numBoundary
    subserialCircle = 0;
    numSegment = geomData(i).numSeg;
    shapeSegment = geomData(i).shape;
    nodeSegment = geomData(i).node;
    circleCenters = geomData(i).circle;
    
    % Count the number of circular segments in the current boundary
    for j = 1:size(shapeSegment, 1)
        if shapeSegment(j) > 0
            subserialCircle = subserialCircle + 1;
        end
    end

    serialSegment = serialSegment + numSegment;

    % Generate NURBS parameters for the current boundary
    nurbsParameters = generateNURBSParameters(numSegment, shapeSegment, ...
        nodeSegment, circleCenters, p(i), isClosed(i));

    % Set refinement level for NURBS parameters
    nurbsParameters.refinement = refinement{i};

    % Generate mesh for the current boundary
    boundaryStr = generateMesh2D(nurbsParameters);

    % Set material properties for the boundary
    boundaryStr.E = E(i);
    boundaryStr.nu = nu;

    % Set number of segments for the boundary
    boundaryStr.numSeg = numSegment;

    % Store boundary structure in array
    nurbsStr(i) = boundaryStr;

    % Plot the mesh
    plotMesh(boundaryStr, numPlot);
end

end
