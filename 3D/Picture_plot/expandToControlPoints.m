function [ctrlDisplacement, ctrlTraction] = expandToControlPoints(nurbsStr, collocNew, disp, trac)
% Expands the displacement and traction values from collocation points to control points.
%
% Inputs:
%   nurbsStr: Structure array containing boundary element data.
%   displacement: Displacement values at control points.
%   traction: Traction values at control points.
%
% Outputs:
%   ctrlDisplacement: Expanded displacement values at control points.
%   ctrlTraction: Expanded traction values at control points.

numPatches = length(nurbsStr);

% Initialize control point displacement and traction
ctrlDisplacement = cell(numPatches, 1);
ctrlTraction = cell(numPatches, 1);
u3 = [disp(1:3:end), disp(2:3:end), disp(3:3:end)];
t3 = [trac(1:3:end), trac(2:3:end), trac(3:3:end)];

[indU, indT, indNumU, indNumT, indNumB] = initializeIndices(nurbsStr, collocNew);

for i = 1:numPatches
    numCtrlPts = size(nurbsStr(i).ctrlPts, 1);
    globU = nurbsStr(i).globU;
    globT = nurbsStr(i).globT;
    b = nurbsStr(i).bouInd;

    % Initialize displacement and traction for this patch
    displacementPatch = zeros(numCtrlPts, 3);
    tractionPatch = zeros(numCtrlPts, 3); % Assuming 3 traction components in 3D

    % Map collocation points to control points
    for j = 1:numCtrlPts
        displacementPatch(j, :) = u3(globU(j) + indNumU(b), :);
    end
    for j = 1:numCtrlPts
        tractionPatch(j, :) = t3(globT(j) + indNumT(b), :);
    end

    ctrlDisplacement{i} = displacementPatch;
    ctrlTraction{i} = tractionPatch;
end

% Combine results from all patches
ctrlDisplacement = cell2mat(ctrlDisplacement);
ctrlTraction = cell2mat(ctrlTraction);
end
