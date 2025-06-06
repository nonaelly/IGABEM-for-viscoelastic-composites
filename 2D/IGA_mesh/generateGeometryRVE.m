% Generate geometry

function geomData = generateGeometryRVE(shapes, params)
% Initialize geometry data structure
geomData = struct('numSeg', [], 'shape', [], 'node', [], 'circle', []);

% Loop over each shape and generate geometry
for i = 1:length(shapes)
    shapeType = shapes{i};
    shapeParams = params{i};
    switch shapeType
        case 'cyl'  % Thick-walled cylinder
            O = shapeParams(1:2);
            R1 = shapeParams(3);
            R2 = shapeParams(4);
            geomData(i) = generateCylinderGeometry(O, R1, R2);
        case 'dumb'  % Dumbbell
            H1 = shapeParams(1);
            H2 = shapeParams(2);
            L1 = shapeParams(3);
            L2 = shapeParams(4);
            L3 = shapeParams(5);
            R = shapeParams(6);
            geomData(i) = generateDumbbellGeometry(H1, H2, L1, L2, L3, R);
        case 'rec'  % rectangle
            O = shapeParams(1:2);
            a = shapeParams(3);
            b = shapeParams(4);
            geomData(i) = generateRectangleGeometry(O, a, b);
        case 'ell'  % Ellipse
            O = shapeParams(1:2);
            a = shapeParams(3);
            if length(shapeParams) == 4
                b = shapeParams(4);
            elseif length(shapeParams) == 3
                b = a;
            else
                error('wrong input for geometry!')
            end
            geomData(i) = generateEllipseGeometry(O, a, b);
            % Add more cases for other shapes as needed
        otherwise
            error('Unsupported shape type');
    end
    if i>1
        % anti-clockwise
        geomData(i).shape = flip(geomData(i).shape);
        geomData(i).node = flip(geomData(i).node);
        geomData(i).circle = flip(geomData(i).circle);
    end
end
end

%% Sub functions.
% =========================================================================
% Function for cylinder
function geomData = generateCylinderGeometry(O, R1, R2)
numSegment = 4;                           % Number of segment.
shapeSegment = [0; 1; 0; 1];              % Shape of segment.
nodeSegment = [R1,0;R2,0;0,R2;0,R1] + O;  % Nodes for each segment
circleCenters = [O; O];                   % Center points for circular arcs.

geomData.numSeg = numSegment;
geomData.shape = shapeSegment;
geomData.node = nodeSegment;
geomData.circle = circleCenters;
end

% Function for Dumbbell
function geomData = generateDumbbellGeometry(H1, H2, L1, L2, L3, R)
numSegment = 6;                           % Number of segment.
shapeSegment = [0; 0; 0; 1; 0; 0];              % Shape of segment.
nodeSegment = [0,0; L1+L2+L3,0; L1+L2+L3,H2; ...
    L1+L2,H2; L1,H1; 0,H1];  % Nodes for each segment
circleCenters = [L1,H1+R];                   % Center points for circular arcs.

geomData.numSeg = numSegment;
geomData.shape = shapeSegment;
geomData.node = nodeSegment;
geomData.circle = circleCenters;
end

% Function for Rectangle
function geomData = generateRectangleGeometry(O, a, b)
numSegment = 4;                           % Number of segment.
shapeSegment = [0; 0; 0; 0];              % Shape of segment.
nodeSegment = [-a,-b;a,-b;a,b;-a,b] + O;  % Nodes for each segment
circleCenters = [];                   % Center points for circular arcs.

geomData.numSeg = numSegment;
geomData.shape = shapeSegment;
geomData.node = nodeSegment;
geomData.circle = circleCenters;
end

% Function for Ellipse
function geomData = generateEllipseGeometry(O, a, b)
numSegment = 4;                           % Number of segment.
shapeSegment = [1; 1; 1; 1];              % Shape of segment.
nodeSegment = [a,0;0,b;-a,0;0,-b] + O;  % Nodes for each segment
circleCenters = [O; O; O; O];                   % Center points for circular arcs.

geomData.numSeg = numSegment;
geomData.shape = shapeSegment;
geomData.node = nodeSegment;
geomData.circle = circleCenters;
end