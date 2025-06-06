function plotMesh(boundaryStructure, numPlot)
% Function: plotMesh
% Description: Plots the NURBS mesh.
%
% Input:
%   boundaryStructure (struct): Struct containing boundary information.
%   numPlot (scalar): Figure number for plotting.

% Extract data from the boundary structure
knotVector = boundaryStructure.knotU;
controlPoints = boundaryStructure.ctrlPts;
collocationPoints = boundaryStructure.collocPts;
weights = boundaryStructure.weights;
p = boundaryStructure.p;
t2uConn = boundaryStructure.t2uConn; 

% Generate parameter values for plotting
n = round(1 / knotVector(p + 2));
xi = linspace(0, max(knotVector), n * 20 + 1);
% xi = linspace(0, max(knotVector), n * 1 + 1);

% Interpolate knot coordinates
knotCoords = zeros(size(xi, 2), 2);
for i = 1:length(knotCoords)
    knotCoords(i, 1) = NURBSinterpolation(xi(i), p, knotVector, controlPoints(:, 1)', weights);
    knotCoords(i, 2) = NURBSinterpolation(xi(i), p, knotVector, controlPoints(:, 2)', weights);
end

% Plot knot coordinates
figure(numPlot);
plot(knotCoords(:, 1), knotCoords(:, 2), '-k', 'LineWidth', 2.5);
hold on

% Plot control points and collocation points
% plot(controlPoints(:, 1), controlPoints(:, 2), 'ro', 'MarkerSize', 14, 'LineWidth', 3)
% plot(collocationPoints(:, 1), collocationPoints(:, 2), 'kx', 'MarkerSize', 14, 'LineWidth', 2.5)
plot(controlPoints(:, 1), controlPoints(:, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2)
plot(collocationPoints(:, 1), collocationPoints(:, 2), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5)

axis equal
title('NURBS Mesh and Collocation Points');
% legend('NURBS','Control Points','Collocation Points', 'Elements', '')
legend('NURBS','Control Points','Collocation Points')
xlabel('X');
ylabel('Y');
hold on
end
