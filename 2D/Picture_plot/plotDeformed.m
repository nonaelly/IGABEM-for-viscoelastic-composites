function plotDeformed(nurbsStr, disp, factor)
% Function: plotDeformed
% Description: Plots the deformed NURBS mesh based on the computed displacements.
%
% Input:
%   nurbsStr (struct): Struct array containing NURBS information.
%   displacement (vector): Computed displacement vector.
%   factor (scalar): Scaling factor for the displacements.

figure;
hold on;

tempCollPts = 0;
for i = 1:length(nurbsStr)
    controlPoints = nurbsStr(i).ctrlPts;
    knotVec = nurbsStr(i).knotU;
    weights = nurbsStr(i).weights;
    p = nurbsStr(i).p;
    
    numCtrlPts = size(controlPoints, 1);
    deformedCtrlPts = zeros(numCtrlPts, 2);

    % Apply displacements to control points
    for j = 1:numCtrlPts-nurbsStr(i).isClosed
        deformedCtrlPts(j, :) = controlPoints(j, 1:2) + ...
            [disp(tempCollPts*2 + j*2-1) disp(tempCollPts*2 + j*2)] * factor;
    end
    tempCollPts = tempCollPts + numCtrlPts - nurbsStr(i).isClosed;
    
    % Ensure closed curves for closed NURBS
    if nurbsStr(i).isClosed
        deformedCtrlPts(numCtrlPts, :) = deformedCtrlPts(1, :);
    end

    % Interpolate and plot the deformed NURBS curve
    numPts = 40*(nurbsStr(i).numSeg)+1;
    xi = linspace(0, max(knotVec), numPts);
    points = zeros(numPts, 2);
    pointsPre = zeros(numPts, 2);
    for point = 1:numPts
        pointsPre(point, 1) = NURBSinterpolation(xi(point), p, knotVec, controlPoints(:, 1)', weights');
        pointsPre(point, 2) = NURBSinterpolation(xi(point), p, knotVec, controlPoints(:, 2)', weights');
        points(point, 1) = NURBSinterpolation(xi(point), p, knotVec, deformedCtrlPts(:, 1)', weights');
        points(point, 2) = NURBSinterpolation(xi(point), p, knotVec, deformedCtrlPts(:, 2)', weights');
    end
    plot(pointsPre(:, 1), pointsPre(:, 2), 'k-');  % Undeformed with solid line
    plot(points(:, 1), points(:, 2), 'r--');  % Deformed with dashed line
end

axis equal;
title('Deformed NURBS Mesh');
xlabel('X');
ylabel('Y');
legend('Undeformed', 'Deformed');
hold off;
end
