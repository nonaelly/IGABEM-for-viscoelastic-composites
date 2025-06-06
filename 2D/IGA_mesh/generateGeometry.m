% =========================================================================
% Generate geometry information for different shapes.
% Supported: 'rec' (rectangle), 'cyl' (cylinder), 'dumb' (dumbbell)
%
% Input:
%   shapes : cell array of shape types (e.g., {'rec','cyl'})
%   params : cell array of parameter vectors corresponding to each shape
%
% Output:
%   geomData : struct array containing numSeg, shape, node, and circle fields
%
% Author: [Wang Zhetong]
% =========================================================================
function geomData=generateGeometry(shapes,params)

% Initialize geometry data structure
geomData=struct('numSeg',[],'shape',[],'node',[],'circle',[]);

% Loop over each shape and generate geometry
for i=1:length(shapes)
    shapeType=shapes{i};
    shapeParams=params{i};
    switch shapeType
        case 'cyl'  % Thick-walled cylinder
            O=shapeParams(1:2);          % Center
            R1=shapeParams(3);           % Inner radius
            R2=shapeParams(4);           % Outer radius
            geomData(i)=generateCylinderGeometry(O,R1,R2);
        case 'rec'  % Rectangle
            O=shapeParams(1:2);          % Bottom-left corner
            a=shapeParams(3);            % Width
            b=shapeParams(4);            % Height
            geomData(i)=generateRectangleGeometry(O,a,b);
        case 'dumb'  % Dumbbell
            H1=shapeParams(1);
            H2=shapeParams(2);
            L1=shapeParams(3);
            L2=shapeParams(4);
            L3=shapeParams(5);
            R=shapeParams(6);
            geomData(i)=generateDumbbellGeometry(H1,H2,L1,L2,L3,R);
        otherwise
            error('Unsupported shape type');
    end
end
end

%% Sub functions.
% =========================================================================
% Function for thick-walled cylinder geometry
function geomData=generateCylinderGeometry(O,R1,R2)
numSegment=4;                               % Number of segments
shapeSegment=[0;1;0;1];                     % Segment types: 1 for arc
nodeSegment=[R1,0;R2,0;0,R2;0,R1]+O;        % Segment nodes
circleCenters=[O;O];                        % Arc centers

geomData.numSeg=numSegment;
geomData.shape=shapeSegment;
geomData.node=nodeSegment;
geomData.circle=circleCenters;
end

% =========================================================================
% Function for rectangular geometry
function geomData=generateRectangleGeometry(O,a,b)
numSegment=4;
shapeSegment=[0;0;0;0];                     % All straight edges
nodeSegment=[0,0;a,0;a,b;0,b]+O;            % Corner coordinates
circleCenters=[O;O];                        % Dummy values for consistency

geomData.numSeg=numSegment;
geomData.shape=shapeSegment;
geomData.node=nodeSegment;
geomData.circle=circleCenters;
end

% =========================================================================
% Function for dumbbell geometry (rectangle with semi-circular notch)
function geomData=generateDumbbellGeometry(H1,H2,L1,L2,L3,R)
numSegment=6;
shapeSegment=[0;0;0;1;0;0];                 % One arc (the notch)
nodeSegment=[0,0;L1+L2+L3,0;L1+L2+L3,H2; ...
              L1+L2,H2;L1,H1;0,H1];         % Nodes defining geometry
circleCenters=[L1,H1+R];                    % Notch arc center

geomData.numSeg=numSegment;
geomData.shape=shapeSegment;
geomData.node=nodeSegment;
geomData.circle=circleCenters;
end

