% =========================================================================
% Generate inner points for specific geometric models
%
% Author: Wang Zhetong
% =========================================================================
function innerPts=getInnerPoints(modelName,varargin)
% Input:
%   modelName : geometry type ("cyl", "2", "dumb")
%   varargin  : {numInner, shapeParameters}
% Output:
%   innerPts  : inner points (n×2) coordinates

numInner=varargin{1};
shapeParameters=varargin{2};

switch modelName
    case "cyl"  % Thick-walled cylinder
        numTheta=numInner(1);           % Number of angles
        numR=numInner(2);               % Number of radii
        O=shapeParameters(1:2);         % Center
        R1=shapeParameters(3);          % Inner radius
        R2=shapeParameters(4);          % Outer radius
        innerPts=zeros(numTheta*numR,2);
        for i=1:numTheta
            for j=1:numR
                R=(R2-R1)*j/(numR+1)+R1;
                theta=90*i/(numTheta+1);  % Angle in degrees
                innerPts((i-1)*numR+j,:)=[R*cosd(theta),R*sind(theta)];
            end
        end
        innerPts=innerPts+O;

    case "dumb"  % Dumbbell-shaped plate
        numX1=numInner(1);              % Region 1 (left)
        numX2=numInner(2);              % Region 2 & 3 (right)
        numY1=numInner(3);              % Height for region 1 & 2
        numY2=numInner(4);              % Height for region 3
        H1=shapeParameters(1);
        H2=shapeParameters(2);
        L1=shapeParameters(3);
        L2=shapeParameters(4);
        L3=shapeParameters(5);
        R=shapeParameters(6);
        totalPoints=numX1*numY1 + numX2*(numY1+1) + numX2*numY2;
        innerPts=zeros(totalPoints,2);

        % Region 1 (rectangle left)
        for i=1:numX1
            for j=1:numY1
                x=i/numX1*L1;
                y=j/(numY1+1)*H1;
                innerPts((i-1)*numY1+j,:)=[x,y];
            end
        end

        % Region 2 (middle rectangle)
        offset=numX1*numY1;
        for i=1:numX2
            for j=1:(numY1+1)
                x=i/(numX2+1)*(L2+L3)+L1;
                y=j/(numY1+1)*H1;
                innerPts(offset+(i-1)*(numY1+1)+j,:)=[x,y];
            end
        end

        % Region 3 (top arc region)
        offset=offset+numX2*(numY1+1);
        for i=1:numX2
            for j=1:numY2
                y=j/(numY2+1)*(H2-H1)+H1;
                xAll=L2+L3-sqrt(R^2-(R-y+H1)^2);
                x=L1+L2+L3-i/(numX2+1)*xAll;
                innerPts(offset+(i-1)*numY2+j,:)=[x,y];
            end
        end

    otherwise
        error('No such case!')
end
end



