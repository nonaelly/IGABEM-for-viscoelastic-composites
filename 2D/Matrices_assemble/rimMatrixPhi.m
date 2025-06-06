function [rimPhi] = rimMatrixPhi(approximatedPoints, appliedPoints, dA, d)
    % Assemble Transform Matrices for RIM.
    %
    % Inputs:
    %   approximatedPoints: The points to be approximated in RIM.
    %   appliedPoints: Applied points in RIM, containing boundary points and inner points.
    %   dA: Support domain size.
    %   d: demension.
    %
    % Outputs:
    %   rimPhi: The transform matrix in radial integration method.
 
    numApplied = size(appliedPoints, 1);
    numApproximated = size(approximatedPoints, 1);
    rimPhi = zeros((numApproximated + d + 1) * d, (numApplied + d + 1) * d);
    for i = 1:numApproximated
        row = i * d - d + 1 : i * d;
        for j = 1:numApplied
            column = j * d - d + 1 : j * d;
            R = norm(approximatedPoints(i, :) - appliedPoints(j, :));
            if R < dA
                phiR = 1 - 6 * ((R / dA)^2) + 8 * ((R / dA)^3) - 3 * ((R / dA)^4);
            else
                phiR = 0;
            end
            [phiSub] = phiR * eye(d);
            rimPhi(row, column) = phiSub;
        end
    end
    
    column = numApplied * d + 1 : (numApplied + 1) * d;
    for i = 1:numApproximated
        row = i * d - d + 1 : i * d;
        [phiSub] = eye(d);
        rimPhi(row, column) = phiSub;
    end
    
    for i = 1:numApproximated    
        row = i * d - d + 1 : i * d;
        for j = numApplied + 2 : numApplied + d + 1
            column = j * d - d + 1 : j * d;
            xK = approximatedPoints(i, j - 1 - numApplied);
            [phiSub] = xK * eye(d);
            rimPhi(row, column) = phiSub;
        end
    end

	row = (numApproximated * d + 1) : (numApproximated + d + 1) * d;
	for i = 1:numApplied
		column = i * d - d + 1 : i * d;
		[phiSub] = [eye(d); eye(d)*appliedPoints(i, 1); eye(d)*appliedPoints(i, 2)];
        if d == 3
            [phiSub] = [phiSub; eye(d)*appliedPoints(i, 3)];
        end
		rimPhi(row, column) = phiSub;
	end
end