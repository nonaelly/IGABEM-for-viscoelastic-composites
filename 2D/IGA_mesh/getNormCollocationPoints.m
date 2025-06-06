function normsOnColloc = getNormCollocationPoints(elemRange, globIdx, ...
    collocXi, UConn, ctrlPts, weights, knotVec, p, isClosed)
% Function: getNormCollocationPoints
% Description: Computes the normal vectors at collocation points on the NURBS curve.
%
% Inputs:
%   elemRange (matrix): Range of elements in terms of knot span indices.
%   globIdx (matrix): Global basis indices of the NURBS basis functions for each element.
%   collocXi (vector): Local parametric coordinates of collocation points.
%   UConn (matrix): Element connectivity for displacement nodes.
%   ctrlPts (matrix): Control points of the NURBS curve.
%   weights (vector): Weights of the control points.
%   knotVec (vector): Knot vector of the NURBS curve.
%   p (scalar): Degree of the NURBS curve.
%   isClosed (logical): Flag indicating whether the NURBS curve is closed.
%
% Output:
%   normsOnColloc (matrix): Normal vectors at collocation points.

% Initialize array to store normal vectors
normsOnColloc = zeros(length(collocXi), 4);

% Define tolerance value
tol = eps;

% Loop over each collocation point
for i = 1:length(collocXi)
    xi = collocXi(i);
    
    % Loop over each element to determine the corresponding element
    for e = 1:size(elemRange, 1)
        range = elemRange(e, :);
        basisIdx = globIdx(e, :);
        elemNodes = UConn(e, :);
        coordElemNodes = ctrlPts(elemNodes, 1:2);
        
        % Check if the collocation point is at the beginning or end of the curve
        if isClosed
            if e == size(elemRange, 1) && xi == 0
                normsOnColloc(i, 1:2) = getNormPoint(coordElemNodes, ...
                    1 - tol, basisIdx, knotVec, weights, p);
            end
        else
            if e == 1 && xi == 0
                normsOnColloc(i, 1:2) = getNormPoint(coordElemNodes, ...
                    xi + tol, basisIdx, knotVec, weights, p);
            elseif e == size(elemRange, 1) && xi == 1
                normsOnColloc(i, 3:4) = getNormPoint(coordElemNodes, ...
                    xi - tol, basisIdx, knotVec, weights, p);
            end
        end
        
        % Check if the collocation point is within the range of the current element
        if xi <= range(2) && xi >= range(1)
            if abs(xi - range(2)) < tol
                normsOnColloc(i, 1:2) = getNormPoint(coordElemNodes, ...
                    xi - tol, basisIdx, knotVec, weights, p);
            elseif abs(xi - range(1)) < tol
                normsOnColloc(i, 3:4) = getNormPoint(coordElemNodes, ...
                    xi, basisIdx, knotVec, weights, p);
            else
                normsOnColloc(i, 1:2) = getNormPoint(coordElemNodes, ...
                    xi, basisIdx, knotVec, weights, p);
                normsOnColloc(i, 3:4) = normsOnColloc(i, 1:2);
            end
        end
    end
end

end
