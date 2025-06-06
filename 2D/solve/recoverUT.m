function [Disp, Trac] = recoverUT(Disp, Trac, U32, U14, T32, T14, HOri, GOri, Q, vecP, PBC, nurbsStr)
% Calculate the u+, t+ and traction on the inclusions' boundaries.
D = 2;
Delta = [0,1; 1,0]';
for i = 1:D
    epsilonDelta = PBC*Delta(:, i);
    epsilonDelta = repmat(epsilonDelta, length(U32{i})/D, 1);
    Disp(U32{i}) = Disp(U14{i}) + epsilonDelta;
    Trac(T32{i}) = -Trac(T14{i});
end

[indU, indT, indNumU, indNumT, indNumB] = initializeIndices2D(nurbsStr);
numU = indNumU(end);
numT = indNumT(end);

% H*uI = G*tI+Q
indInT = nurbsStr(1).numCollocT*D+1 : numT*D;
indInU = nurbsStr(1).numCollocU*D+1 : numU*D;
indOutT = 1:nurbsStr(1).numCollocT*D;

Trac(indInT) = GOri(indInU, indInT)\(HOri(indInU,:)*Disp - ...
    GOri(indInU, indOutT)*Trac(indOutT) - Q(indInU) + vecP(indInU));
end