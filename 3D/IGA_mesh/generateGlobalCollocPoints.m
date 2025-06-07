function [nurbsStr,collocNew] = generateGlobalCollocPoints(nurbsStr,numPlot)
    % Generates global collocation points for all boundaries and records global indices.
    %
    % Input:
    %   nurbsStr (struct array): Array of structs containing NURBS data for each patch.
    %
    % Output:
    %   collocNew (struct): Struct containing new collocation points and global indices.

    % Initialize variables
    tol = 1e-5;
    collocNew = [];
    indPatch = 1;
    tempT = 0;

    % Generate new collocation points for each boundary
    for b = 1:max([nurbsStr.bouInd])
        boundaryPatches = nurbsStr([nurbsStr.bouInd] == b);
        collocPtsNew = [];
        infoPts = [];
        for j = 1:length(boundaryPatches)
            collocPts = boundaryPatches(j).collocPts;
            numCollocPts = size(collocPts, 1);
            locToGloIndex = zeros(1, numCollocPts);
            
            for k = 1:numCollocPts
                temp_c = collocPts(k, :);
                % Check if the point already exists in collocPtsNew
                [isRepeated, repeatedIndex] = isRepeatedCoord(temp_c, collocPtsNew, tol);
                if ~isRepeated
                    collocPtsNew = [collocPtsNew; temp_c];
                    locToGloIndex(k) = size(collocPtsNew, 1);
                    infoPts = [infoPts; j,k];
                else
                    locToGloIndex(k) = repeatedIndex;
                end
            end
            boundaryPatches(j).locToGloIndex = locToGloIndex;
            nurbsStr(indPatch).globU = locToGloIndex;
            if b == 1 % Cube
                nurbsStr(indPatch).globT = 1+tempT :boundaryPatches(j).numCollocT +tempT;
                tempT = tempT + boundaryPatches(j).numCollocT;
            else
                nurbsStr(indPatch).globT = locToGloIndex;
            end  
            indPatch = indPatch + 1;
        end
        collocNew(b).collocPts = collocPtsNew;
        collocNew(b).infoPts = infoPts;
        collocNew(b).boundaryIndex = b;

        % Plot the collocation points for each boundary
        figure(numPlot);
        hold on;
        if b>1
%         for i = 1:length(collocPtsNew)
%             plot3(collocPtsNew(i, 1), collocPtsNew(i, 2), collocPtsNew(i, 3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'b');
%             text(collocPtsNew(i, 1)-0.01, collocPtsNew(i, 2)-0.01, collocPtsNew(i, 3), num2str(i), 'FontSize', 8, ...
%                 'Color', 'k', 'HorizontalAlignment', 'right');
%         end
        end
    end
    
    axis equal;
    hold on;
end

%======================================================
function [isRepeated, repeatedIndex] = isRepeatedCoord(coord, coordList, tol)
    % Checks if a coordinate already exists in a list of coordinates within a tolerance.
    %
    % Input:
    %   coord (vector): Coordinate to check.
    %   coordList (matrix): List of existing coordinates.
    %   tol (scalar): Tolerance for coordinate comparison.
    %
    % Output:
    %   isRepeated (logical): True if the coordinate is repeated, false otherwise.
    %   repeatedIndex (scalar): Index of the repeated coordinate in the list.

    isRepeated = false;
    repeatedIndex = 0;
    for i = 1:size(coordList, 1)
        if norm(coord - coordList(i, :)) < tol
            isRepeated = true;
            repeatedIndex = i;
            break;
        end
    end
end
