function [H, G, vecP] = AssembleMatrix(HElaBou, GElaBou, tranU, tranT, U32, U14, T32, T14, PBC)

D = 2;

H = HElaBou/tranU; 
G = GElaBou/tranT;

lenU3 = length(U32{1});
lenU2 = length(U32{2});
% u3 - u1
delta31 = PBC*[0;1];
% u2 - u4
delta24 = PBC*[1;0];
deltaX = [repmat(delta31, lenU3/2, 1); repmat(delta24, lenU2/2, 1)];
vecP = - H(:, [U32{1}; U32{2}])*deltaX;

% Rebulid the Matrix G and H.
H(:, [U14{1}; U14{2}]) = H(:, [U14{1}; U14{2}]) + H(:, [U32{1}; U32{2}]);
G(:, [T14{1}; T14{2}]) = G(:, [T14{1}; T14{2}]) - G(:, [T32{1}; T32{2}]);

end