function [H, G, vecP] = rebulidPeriodicMatrix(H, G, face325U, face146U, face325T, face146T, epsilonPB)

D = 3;
deltaXs = {[0; 1; 0], [1; 0; 0], [0; 0; 1]};

% {u+} = {u-} + [ε]*{Δx}; {t+} = -{t-} 
% [H+, H-]*{u+; u-} = [G+, G-]*{t+; t-} + {Q} 

vecP = zeros(size(H,1),1);
for i = 1:D
    n = size(face146U{i},1)/3;
    deltaX = repmat(epsilonPB*deltaXs{i}, n, 1);
    H(:, face146U{i}) = H(:, face146U{i}) + H(:, face325U{i});
    G(:, face146T{i}) = G(:, face146T{i}) - G(:, face325T{i});
    vecP = - H(:, face325U{i})*deltaX + vecP;
end


end