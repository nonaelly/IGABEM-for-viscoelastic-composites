function plotMesh3D(boundaryStr, numPlot)
% plotMesh3D Plots the 3D NURBS mesh.
%
% Inputs:
%   boundaryStructure (struct): Struct containing boundary information.
%   numPlot (scalar): Figure number for plotting.

% Create figure for plotting
figure(numPlot)
hold on

% Loop over each boundary structure
for i = 1:length(boundaryStr)
    ctrlPts = boundaryStr(i).ctrlPts;
    knotU = boundaryStr(i).knotU;
    knotV = boundaryStr(i).knotV;
    collocPts = boundaryStr(i).collocPts;
    
    % Plot control points
    plot3(ctrlPts(:, 1), ctrlPts(:, 2), ctrlPts(:, 3), 'bo', 'MarkerSize', 3);
    
    % Plot control point mesh in U direction
    numCtrlPtsU = length(knotU) - boundaryStr(i).p - 1;
    numCtrlPtsV = length(knotV) - boundaryStr(i).q - 1;
    
    for k = 1:numCtrlPtsV
        ind = (1:numCtrlPtsU) + (k - 1) * numCtrlPtsU;
        plot3(ctrlPts(ind, 1), ctrlPts(ind, 2), ctrlPts(ind, 3), 'b--');
    end
    
    % Plot control point mesh in V direction
    for j = 1:numCtrlPtsU
        ind = (1:numCtrlPtsU:numCtrlPtsU * numCtrlPtsV) + j - 1;
        plot3(ctrlPts(ind, 1), ctrlPts(ind, 2), ctrlPts(ind, 3), 'b--');
    end
    
    % Plot collocation points
    plot3(collocPts(:, 1), collocPts(:, 2), collocPts(:, 3), 'ro', 'MarkerSize', 6, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'b');
    
    % Optionally label collocation points
%     if boundaryStr(i).refinement < 4
%         for j = 1:size(collocPts, 1)
%             text(collocPts(j, 1)-0.01, collocPts(j, 2)-0.01, collocPts(j, 3), num2str(j), 'FontSize', 8, ...
%                 'Color', 'k', 'HorizontalAlignment', 'right');
%         end
%     end
end

% Set plot attributes
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
title('3D NURBS Mesh')
grid on
view(-37.5 + 90, 30)
% hold off
end
