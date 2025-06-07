function controlPoints = getPointsCube(lengthCube, centerCube, patch)
% Generates control points for a cubic patch
% 
% Inputs:
% - lengthCube: Vector [a, b, c] defining the dimensions of the cube
% - centerCube: Vector [x, y, z] defining the center of the cube
% - patch: Specifies which patch of the cube to generate points for (1 to 6)
% 
% Output:
% - controlPoints: 9x3 matrix of control points
% 
%           ______
%          /     /| →3
%         /  5  / |
%        /_____/ 2|
%     4← |     |  /
%        |  1  | /      z  
%        |_____|/       ↑   y
%           ↑           |  ↗ 
%           6           O---→ x

% -------------------------------------------------------------------------
a = lengthCube(1); 
b = lengthCube(2); 
c = lengthCube(3);

controlPoints = zeros(9, 3);

switch patch
    case {1, 3}
        if patch == 1
            controlPoints(:, 1) = a * [-1; 0; 1; -1; 0; 1; -1; 0; 1];
            controlPoints(:, 2) = -b * ones(9, 1);
        else
            controlPoints(:, 1) = -a * [-1; 0; 1; -1; 0; 1; -1; 0; 1];
            controlPoints(:, 2) = b * ones(9, 1);
        end
        controlPoints(:, 3) = [-c; -c; -c; 0; 0; 0; c; c; c];
        
    case {2, 4}
        if patch == 4
            controlPoints(:, 2) = -b * [-1; 0; 1; -1; 0; 1; -1; 0; 1];
            controlPoints(:, 1) = -a * ones(9, 1);
        else
            controlPoints(:, 2) = b * [-1; 0; 1; -1; 0; 1; -1; 0; 1];
            controlPoints(:, 1) = a * ones(9, 1);
        end
        controlPoints(:, 3) = [-c; -c; -c; 0; 0; 0; c; c; c];

    case {5, 6}
        controlPoints(:, 1) = a * [-1; 0; 1; -1; 0; 1; -1; 0; 1];
        if patch == 5
            controlPoints(:, 2) = b * [-1; -1; -1; 0; 0; 0; 1; 1; 1];
            controlPoints(:, 3) = c * ones(9, 1);
        else
            controlPoints(:, 2) = -b * [-1; -1; -1; 0; 0; 0; 1; 1; 1];
            controlPoints(:, 3) = -c * ones(9, 1);
        end

    otherwise
       error('Wrong Patch')
end

% Translate the control points to the center of the cube
for j = 1:3
    controlPoints(:, j) = controlPoints(:, j) + centerCube(j);
end

end
