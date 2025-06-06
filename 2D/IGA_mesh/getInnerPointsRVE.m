function innerPoints = getInnerPointsRVE(modelName, varargin)
% Generates inner points within a Representative Volume Element (RVE) for composite materials.
%
% Input:
%   modelName : Indicates the type of model. For RVE, only "compositeRVE" is valid.
%   varargin : Additional variable-length input arguments.
%       numInner (array): Number of inner points in each direction for the rectangular region.
%       shapeParameters (cell array): Parameters defining the composite RVE:
%           numInnerCircle (array): Number of inner points in radial and angular directions for the circular regions.
%           circleCenters (matrix): Coordinates of the centers of the circles.
%           circleRadii (matrix): Radii of the circles or ellipses.
%           lengthMatrix (array): Lengths of the RVE in x and y directions.
%           centerMatrix (array): Coordinates of the center of the RVE.
%           scale (scalar): Scaling factor for the circles/ellipses.
%           tolerance (scalar): Tolerance for determining whether a point is inside a circle/ellipse.
%
% Output:
%   innerPoints (matrix): Matrix containing the coordinates of the inner points.

% Extract input arguments
numInner = varargin{1};
shapeParameters = varargin{2};

switch modelName
    case "RVE" % Only case for RVE
        % Extract parameters for RVE
        numInnerX = numInner(1);
        numInnerY = numInner(2);
        numR = numInner(3);
        numAgl = numInner(4);
        scale = numInner(5);
        tolerance = numInner(6);

        numInc = length(shapeParameters) - 1;
        [circleCenters, circleRadii] = deal(zeros(numInc,2));
        for i = 1:numInc
            circleCenters(i,:) = shapeParameters{i+1}(1:2);
            circleRadii(i,:) = shapeParameters{i+1}(3:4);
        end
        lengthMatrix = shapeParameters{1}(3:4);
        centerMatrix = shapeParameters{1}(1:2);


        % Initialize the matrix for storing inner points
        totalPoints = 2 * numInnerX * numInnerY + 2 * numR * numAgl * size(circleRadii, 1);
        innerPoints = zeros(totalPoints, 2);
        count = 1;
        ox = centerMatrix(1);
        oy = centerMatrix(2);
        a = lengthMatrix(1);
        b = lengthMatrix(2);

        % Generate points within the rectangular region of the RVE
        for ix = 1:numInnerX
            for iy = 1:numInnerY
                x = (ix - 0.5) * 2 * a / numInnerX - a + ox;
                y = (iy - 0.5) * 2 * b / numInnerY - b + oy;
                % Check if the point is inside any circle/ellipse
                if ~isInsideAnyCircle(x, y, circleCenters, circleRadii * (1 + scale), tolerance)
                    innerPoints(count, :) = [x, y];
                    count = count + 1;
                end
            end
        end

        % Generate points within the circular regions of the RVE
        for i = 1:size(circleRadii, 1)
            for ix = 1:numR
                for iy = 1:numAgl
                    x = cos(2 * pi / numAgl * iy) * (ix / numR * circleRadii(i, 1) * scale + circleRadii(i, 1)) + circleCenters(i, 1);
                    y = sin(2 * pi / numAgl * iy) * (ix / numR * circleRadii(i, 2) * scale + circleRadii(i, 2)) + circleCenters(i, 2);
                    innerPoints(count, :) = [x, y];
                    count = count + 1;
                end
            end
        end

        % Remove any unused rows in the innerPoints matrix
        innerPoints = innerPoints(1:count-1, :);

    otherwise
        error('Invalid modelName. Use "RVE" for RVE generation.')
end

end

function inside = isInsideAnyCircle(x, y, circleCenters, circleRadii, tolerance)
% Checks if the point (x, y) is inside any of the circles/ellipses.
%
% Input:
%   x, y : Coordinates of the point to check.
%   circleCenters (matrix): Coordinates of the centers of the circles/ellipses.
%   circleRadii (matrix): Radii of the circles/ellipses.
%   tolerance (scalar): Tolerance for determining proximity to the boundary.
%
% Output:
%   inside (logical): True if the point is inside any circle/ellipse, false otherwise.

inside = false;
for i = 1:size(circleCenters, 1)
    dx = x - circleCenters(i, 1);
    dy = y - circleCenters(i, 2);
    % Determine if the point is inside the ellipse or circle
    if size(circleRadii, 2) == 2
        distance = (dx / circleRadii(i, 1))^2 + (dy / circleRadii(i, 2))^2 - tolerance;
    else
        distance = (dx / circleRadii(i))^2 + (dy / circleRadii(i))^2 - tolerance;
    end

    if distance < 1
        inside = true;
        break;
    end
end

end
