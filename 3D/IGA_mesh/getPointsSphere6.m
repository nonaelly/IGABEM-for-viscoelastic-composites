function [controlPoints, weights, knotU, knotV] = getPointsSphere6(r, o)
% Generates control points, weights, and knot vectors for a six-patch NURBS sphere.
%
% Input:
%   r (scalar): Radius of the sphere.
%   o (vector): Center of the sphere [x, y, z].
%
% Output:
%   controlPoints (cell array): Control points for each patch.
%   weights (cell array): Weights for each patch.
%   knotU (vector): Knot vector in U direction.
%   knotV (vector): Knot vector in V direction.

% Define control points and weights for each of the six patches
controlPoints = cell(6, 1);
weights = cell(6, 1);

% Knot vectors for a 4th degree NURBS
knotU = [0 0 0 0 0 1 1 1 1 1];
knotV = [0 0 0 0 0 1 1 1 1 1];

% Helper variables
sq2 = sqrt(2);
sq3 = sqrt(3);

% Define control points and weights for the first patch
% Initial patch
controlPoints{6} = [
    4*(1-sq3), -4*(1-sq3), 4*(1-sq3);
    -sq2, -sq2*(sq3-4), sq2*(sq3-4);
    0, -4*(1-2*sq3)/3, 4*(1-2*sq3)/3;
    sq2, -sq2*(sq3-4), sq2*(sq3-4);
    -4*(1-sq3), -4*(1-sq3), 4*(1-sq3);

    sq2*(sq3-4), sq2, sq2*(sq3-4);
    (2-3*sq3)/2, -(2-3*sq3)/2, -(sq3+6)/2;
    0, -sq2*(2*sq3-7)/3, -5/3*sq2*sq3;
    -(2-3*sq3)/2, -(2-3*sq3)/2, -(sq3+6)/2;
    -sq2*(sq3-4), sq2, sq2*(sq3-4);

    4*(1-2*sq3)/3, 0, 4*(1-2*sq3)/3;
    sq2*(2*sq3-7)/3, 0, -5*sq2*sq3/3;
    0, 0, 4*(sq3-5)/3;
    -sq2*(2*sq3-7)/3, 0, -5*sq2*sq3/3;
    -4*(1-2*sq3)/3, 0, 4*(1-2*sq3)/3;

    sq2*(sq3-4), -sq2, sq2*(sq3-4);
    (2-3*sq3)/2, (2-3*sq3)/2, -(sq3+6)/2;
    0, sq2*(2*sq3-7)/3, -5/3*sq2*sq3;
    -(2-3*sq3)/2, (2-3*sq3)/2, -(sq3+6)/2;
    -sq2*(sq3-4), -sq2, sq2*(sq3-4);

    4*(1-sq3), 4*(1-sq3), 4*(1-sq3);
    -sq2, sq2*(sq3-4), sq2*(sq3-4);
    0, 4*(1-2*sq3)/3, 4*(1-2*sq3)/3;
    sq2, sq2*(sq3-4), sq2*(sq3-4);
    -4*(1-sq3), 4*(1-sq3), 4*(1-sq3);
    ];

weights{6} = [
    4*(3-sq3); sq2*(3*sq3-2); 4*(5-sq3)/3; sq2*(3*sq3-2); 4*(3-sq3);
    sq2*(3*sq3-2); (sq3+6)/2; sq2*(sq3+6)/3; (sq3+6)/2; sq2*(3*sq3-2);
    4*(5-sq3)/3; sq2*(sq3+6)/3; 4*(5*sq3-1)/9; sq2*(sq3+6)/3; 4*(5-sq3)/3;
    sq2*(3*sq3-2); (sq3+6)/2; sq2*(sq3+6)/3; (sq3+6)/2; sq2*(3*sq3-2);
    4*(3-sq3); sq2*(3*sq3-2); 4*(5-sq3)/3; sq2*(3*sq3-2); 4*(3-sq3)
    ];

controlPoints{6} = controlPoints{6}./weights{6};

% Rotate the initial patch to create the other five patches
rotationMatrices = {
    rot('x',-90), rot('z',90), rot('z',180), rot('z',-90), rot('x',180)
    };

for i = [1,5]
    controlPoints{i} = (rotationMatrices{i} * controlPoints{6}')';
    weights{i} = weights{6};
end

for i = 2:4
    controlPoints{i} = (rotationMatrices{i} * controlPoints{1}')';
    weights{i} = weights{6};
end

for i = 1:6
    controlPoints{i} = r * controlPoints{i} + o;
end

end

%======================================================
function R = rot(axis,theta)
% Rotation matrix for rotation around axis by theta degrees (anti-clock)
switch axis
    case 'x'
        R = [
            1 0 0;
            0 cosd(theta) -sind(theta);
            0 sind(theta) cosd(theta)
            ];
    case 'y'
        R = [
            cosd(theta) 0 sind(theta);
            0 1 0;
            -sind(theta) 0 cosd(theta)
            ];
    case 'z'
        R = [
            cosd(theta) -sind(theta) 0;
            sind(theta) cosd(theta) 0;
            0 0 1
            ];
    otherwise
        error('no such axis!')
end

end

