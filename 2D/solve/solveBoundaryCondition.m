function [displacement, traction] = solveBoundaryCondition(H, G, Q, ...
    knownU, knownT, unknownU, unknownT, collocT, collocU)
% Solves the boundary conditions using the Hu = Gt + Q formulation.

% Inputs:
% - H: Matrix for displacements.
% - G: Matrix for tractions.
% - Q: External load vector.
% - knownU: Indices of known displacements.
% - knownT: Indices of known tractions.
% - unknownU: Indices of unknown displacements.
% - unknownT: Indices of unknown tractions.
% - collocT: Traction values at collocation points.
% - collocU: Displacement values at collocation points.

% Outputs:
% - displacement: Computed displacement vector.
% - traction: Computed traction vector.

% -------------------------------------------------------------------------
% Assemble the matrix.

% Extract known and unknown parts of the matrices.
HKnownU = H(:, knownU);
HUnknownU = H(:, unknownU);
GUnknownT = G(:, unknownT);
GKnownT = G(:, knownT);

% Construct the combined matrices for solving.
A = [-GUnknownT, HUnknownU];
W = [-HKnownU, GKnownT];

% Known traction and displacement values.
knownTraction = collocT(knownT); % Known traction.
knownDisplacement = collocU(knownU); % Known displacement.

% Compute the right-hand side.
z = W * [knownDisplacement; knownTraction] + Q;

% -------------------------------------------------------------------------
% Solve
if size(A, 1) ~= size(A, 2)
    warning('Matrix A is not square. Using pseudoinverse to solve.');
    soln = pinv(A) * z;
else
    % Solve the system of equations.
    soln = A \ z;
end

% Initialize the traction and displacement vectors with collocation values.
traction = collocT;
displacement = collocU;

% Update unknown traction and displacement values.
traction(unknownT) = soln(1:length(unknownT));
displacement(unknownU) = soln((length(unknownT) + 1):end);

end
