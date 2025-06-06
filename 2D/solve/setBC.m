function [preU, preT, unU, unT, vecT, vecU] = setBC(BC,nurbsStr,pressure,dirichlet)
% Generate boundary condition for each 2D boundary.
D = 2;
numU = nurbsStr.numCollocU;
numT = nurbsStr.numCollocT;
vecU = zeros(D*numU,1);
vecT = zeros(D*numT,1);
preU = [];
preT = [];
indDir = 0;
indPre = 0;

if nurbsStr.numSeg ~= length(BC)
    error('BC condition does not match')
end

temp = 0;
for i=1:length(BC)
    uConn = nurbsStr.uConn;
    tConn = nurbsStr.tConn;
    % Index of elements on the segment.
    n = nurbsStr.refine(i);
    indElem = 1+temp : n+temp;
    temp = temp + n;
    switch BC(i)
        case {'x','dx'}
            tempU = uConn(indElem,:)*D-1;
            preU = [preU; reshape(tempU,[],1)];
            tempT = tConn(indElem,:)*D;
            preT = [preT; reshape(tempT,[],1)];
            if strcmpi(BC(i),'x') 
                indDir = indDir + 1;
                vecU(tempU) = dirichlet(indDir);
            end
        case {'y','dy'}
            tempU = uConn(indElem,:)*D;
            preU = [preU; reshape(tempU,[],1)];
            tempT = tConn(indElem,:)*D-1;
            preT = [preT; reshape(tempT,[],1)];
            if strcmpi(BC(i),'y') 
                indDir = indDir + 1;
                vecU(tempU) = dirichlet(indDir);
            end
        case {'f','p'}
            tempT = [tConn(indElem,:)*D-1; tConn(indElem,:)*D];
            preT = [preT; reshape(tempT,[],1)];
            if strcmpi(BC(i),'p') 
                indPre = indPre + 1;
                tempP = unique(reshape(tConn(indElem,:),1,[]));
                normP = nurbsStr.normCollocT(tempP,:);
                vecT(tempP*D-1) = -pressure(indPre) * normP(:,1);
                vecT(tempP*D) = -pressure(indPre) * normP(:,2);
            end
        case 'd' % Fix
            tempU = [uConn(indElem,:)*D-1; uConn(indElem,:)*D];
            preU = [preU; reshape(tempU,[],1)];
        otherwise
            error('Wrong BC')
    end
end

preU = unique(preU,'sorted');
preT = unique(preT,'sorted');
unU = setdiff((1:D*numU)',preU);
unT = setdiff((1:D*numT)',preT);
end