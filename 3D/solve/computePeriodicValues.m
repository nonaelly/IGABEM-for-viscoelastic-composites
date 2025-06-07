function [disp, trac] = computePeriodicValues(disp, trac, face325U, face146U, ...
    face325T, face146T, epsilonPB, nurbsStr, collocNew, HO, GO, Q, vecP)
% {u+} = {u-} + [ε]*{Δx}; {t+} = -{t-}; 
D = 3;
deltaXs = {[0; 1; 0], [1; 0; 0], [0; 0; 1]};

% Initialize
[indU, indT, indNumU, indNumT, ~] = initializeIndices(nurbsStr, collocNew);
numU = indNumU(end);
numT = indNumT(end);

for i = 1:D
    n = size(face146U{i},1)/3;
    deltaX = repmat(epsilonPB*deltaXs{i}, n, 1);
    disp(face325U{i}) = disp(face146U{i}) + deltaX;
    trac(face325T{i}) = - trac(face146T{i});
end

numRVET = sum(indT(:,1) == 1);
indIncluT = 1+numRVET*D : numT*D;
indRVET = 1 : D*numRVET;
numRVEU = sum(indU(:,1) == 1);
indIncluU = 1+numRVEU*D : numU*D;

% [HO]*{u} - [GO(out)]*{t(out)} - {Q} + {Periodic} =  [GO(in)]*{t(in)}
% traction on inclusions' boundaries.
vecLeft = (HO(indIncluU,:)*disp - GO(indIncluU,indRVET)*trac(indRVET) ...
    - Q(indIncluU,:) + vecP(indIncluU,:));
trac(indIncluT) = GO(indIncluU,indIncluT)\ vecLeft;

end