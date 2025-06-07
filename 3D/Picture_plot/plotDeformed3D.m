function plotDeformed3D(nurbsStr, disp, factor, numPlot)
% Function: plotDeformed3D
% Description: Plots the deformed NURBS mesh in 3D based on the computed displacements.
%
% Input:
%   nurbsStr (struct): Struct array containing NURBS information.
%   disp (vector): Computed displacement vector.
%   factor (scalar): Scaling factor for the displacements.
%   numPlot (scalar): Figure number.

figure(numPlot);
hold on;

tempCollPts = 0;
for i = 1:length(nurbsStr)
    controlPoints = nurbsStr(i).ctrlPts; % Now includes [x, y, z]
    knotVecU = nurbsStr(i).knotU;
    knotVecV = nurbsStr(i).knotV; % New for 3D version
    weights = nurbsStr(i).weights;
    p = nurbsStr(i).p;
    q = nurbsStr(i).q; % New for 3D version
    neU = length(unique(knotVecU)) - 1;
    neV = length(unique(knotVecV)) - 1;
    
    numCtrlPts = size(controlPoints, 1);
    deformedCtrlPts = zeros(numCtrlPts, 3);

    % Apply displacements to control points (3D version)
    for j = 1:numCtrlPts
        deformedCtrlPts(j, :) = controlPoints(j, 1:3) + disp(tempCollPts + j, 1:3) * factor;
    end
    tempCollPts = tempCollPts + numCtrlPts;

    % Interpolate and plot the deformed NURBS surface
    num = 10;
    numPtsU = num * neU + 1;
    numPtsV = num * neV + 1;
    xiU = linspace(0, max(knotVecU), numPtsU);
    xiV = linspace(0, max(knotVecV), numPtsV);

    points = zeros(numPtsU, numPtsV, 3);
    pointsPre = zeros(numPtsU, numPtsV, 3);
    for u = 1:numPtsU
        for v = 1:numPtsV
            pointsPre(u, v, :) = NURBSinterpolation3d(xiU(u), xiV(v), p, q, knotVecU, knotVecV, controlPoints, weights');
            points(u, v, :) = NURBSinterpolation3d(xiU(u), xiV(v), p, q, knotVecU, knotVecV, deformedCtrlPts, weights');
        end
    end

    % Plot only the boundary lines of the undeformed and deformed surfaces
    % Plot edges along xiU direction
    for v = [1, numPtsV]  % Only plot the first and last rows for boundary edges
        plot3(pointsPre(:, v, 1), pointsPre(:, v, 2), pointsPre(:, v, 3), 'k-'); % Undeformed in solid black
        plot3(points(:, v, 1), points(:, v, 2), points(:, v, 3), 'r--'); % Deformed in dashed red
    end
    % Plot edges along xiV direction
    for u = [1, numPtsU]  % Only plot the first and last columns for boundary edges
        plot3(pointsPre(u, :, 1), pointsPre(u, :, 2), pointsPre(u, :, 3), 'k-'); % Undeformed in solid black
        plot3(points(u, :, 1), points(u, :, 2), points(u, :, 3), 'r--'); % Deformed in dashed red
    end
end

axis equal
view(-37.5 + 90, 30)
title('Deformed NURBS Surface');
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('Undeformed', 'Deformed');
hold off;
end
