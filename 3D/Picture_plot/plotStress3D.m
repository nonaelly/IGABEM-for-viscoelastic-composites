function plotStress3D(nurbsStr, ctrlStress, stressComponent, figureNum, range)
% Plots the specified stress component on the 3D boundary structures.
%
% Inputs:
%   nurbsStr: Structure array containing boundary element data.
%   ctrlStress: Control point stresses (n*3 vector).
%   stressComponent: The stress component to plot ('xx', 'yy', 'zz', 'xy', 'yz', 'zx').
%   figureNum: Figure number for the plot.

figure(figureNum);
hold on;

numPatches = length(nurbsStr);
indNumT = zeros(1 + numPatches, 1);
for i = 1:numPatches
    indNumT(i + 1) = indNumT(i) + nurbsStr(i).numCollocT;
end

for patchOrder = range
    pointIndex = 1 + indNumT(patchOrder) : indNumT(patchOrder + 1);
    controlPts = nurbsStr(patchOrder).ctrlPts;
    weights = nurbsStr(patchOrder).weights;
    knotVecU = nurbsStr(patchOrder).knotU;
    knotVecV = nurbsStr(patchOrder).knotV;
    p = nurbsStr(patchOrder).p;
    q = nurbsStr(patchOrder).q;
    neU = length(unique(knotVecU)) - 1;
    neV = length(unique(knotVecV)) - 1;

    stress_patch = ctrlStress(pointIndex, :);
    % stress_patch = ctrlStress(nurbsStr(patchOrder).globU, :); you can use
    % this when you input sA as ctrlStress.

    num = 20;
    numU = neU * num;
    numV = neV * num;
    xiU = linspace(0, max(knotVecU), numU + 1);
    xiV = linspace(0, max(knotVecV), numV + 1);

    coords = zeros((numU + 1) * (numV + 1), 3);
    for row = 1:numV + 1
        for col = 1:numU + 1
            coords(col + (row - 1) * (numU + 1), :) = NURBSinterpolation3d(xiU(col), xiV(row), p, q, knotVecU, knotVecV, controlPts, weights');
        end
    end

    stressVals = zeros((numU + 1) * (numV + 1), 1);
    for row = 1:numV + 1
        for col = 1:numU + 1
            stressAtPoint = zeros(6,1);
            stressAtPoint(1:3) = NURBSinterpolation3d(xiU(col), xiV(row), ...
                p, q, knotVecU, knotVecV, stress_patch(:,1:3), weights');
            stressAtPoint(4:6) = NURBSinterpolation3d(xiU(col), xiV(row), ...
                p, q, knotVecU, knotVecV, stress_patch(:,4:6), weights');
            switch stressComponent
                case 'xx'
                    stressVals(col + (row - 1) * (numU + 1)) = stressAtPoint(1);
                case 'yy'
                    stressVals(col + (row - 1) * (numU + 1)) = stressAtPoint(2);
                case 'zz'
                    stressVals(col + (row - 1) * (numU + 1)) = stressAtPoint(3);
                case 'xy'
                    stressVals(col + (row - 1) * (numU + 1)) = stressAtPoint(4);
                case 'yz'
                    stressVals(col + (row - 1) * (numU + 1)) = stressAtPoint(5);
                case 'zx'
                    stressVals(col + (row - 1) * (numU + 1)) = stressAtPoint(6);
                otherwise
                    error('Invalid stress component specified');
            end
        end
    end

    f = zeros(numU * numV, 4);
    for row = 1:numV
        for col = 1:numU
            f(col + (row - 1) * numU, :) = [col col + numU + 1 col + 1 + numU + 1 col + 1] + (row - 1) * (numU + 1);
        end
    end
    S.Faces = f;
    S.Vertices = coords;
    S.FaceVertexCData = stressVals;
    S.FaceColor = 'interp';
    S.EdgeColor = 'none';
    patch(S);

    if patchOrder == range(1)
        maxVal = max(stressVals);
        minVal = min(stressVals);
    else
        maxVal = max(max(stressVals), maxVal);
        minVal = min(min(stressVals), minVal);
    end
end
xlabel('X')
ylabel('Y')
zlabel('Z')
title({['Max value: ', num2str(maxVal)], ['Min value: ', num2str(minVal)]})
colormap jet
colorbar
axis equal
view(-37.5 + 90, 30)
end
