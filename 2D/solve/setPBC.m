function [preU, preT, unU, unT, vecT, vecU, U32, U14, T32, T14] = setPBC(PBC,nurbsStr)
% Periodical boundary condition for each 2D boundary.
% Periodical boundary; 1-3, 4-2.
%      3
%   D-<<-C
% 4 |    | 2
%   |    |
%   A->>-B
%      1

dem = 2;
% Initialize
[indU, indT, indNumU, indNumT, indNumB] = initializeIndices2D(nurbsStr);
numU = indNumU(end);
numT = indNumT(end);

vecU = zeros(dem*numU,1);
vecT = zeros(dem*numT,1);
preU = [];
preT = [];
[U32, U14, T32, T14] = deal(cell(2,1));

for i=1:length(nurbsStr)
    temp = 0;
    for j = 1:nurbsStr(i).numSeg
        uConn = nurbsStr(i).uConn;
        tConn = nurbsStr(i).tConn;
        % Index of elements on the segment.
        n = nurbsStr(i).refine(j);
        indElem = 1+temp : n+temp;
        temp = temp + n;
        if i == 1 % Outer boundary
            tempU = dem*indNumU(i) + [uConn(indElem,:)*dem-1; uConn(indElem,:)*dem];
            tempT = dem*indNumT(i) + [tConn(indElem,:)*dem-1; tConn(indElem,:)*dem];
            tempUP = reshape(uConn(indElem,:),[],1);
            tempUP = tempUP(2:end-1);
            tempTP = reshape(tConn(indElem,:),[],1);
            switch j
                case {1,4}
                    if j==1
                        U14{1} = tempUP;
                        T14{1} = tempTP;
                    else
                        U14{2} = tempUP;
                        T14{2} = tempTP;
                    end
                case {2,3}
                    if j==3
                        U32{1} = flip(tempUP);
                        T32{1} = flip(tempTP);
                    else
                        U32{2} = flip(tempUP);
                        T32{2} = flip(tempTP);
                    end
                    preU = [preU; reshape(tempU,[],1)];
                    preT = [preT; reshape(tempT,[],1)];
            end
        else
            tempT = dem*indNumT(i) + [tConn(indElem,:)*dem-1; tConn(indElem,:)*dem];
            preT = [preT; reshape(tempT,[],1)];
        end
    end

    if i == 1
        % We need to consider the corner points
        % fix (0,0) to avoid rigid displacement.
        A = find(vecnorm(nurbsStr(i).collocPts-[0,0], Inf, 2)<1e-9);
        if size(A)>1
            error('More than one (0,0) !')
        end
        preU = [preU; A*2-1; A*2];

        B = find(vecnorm(nurbsStr(i).collocPts-[1,0], Inf, 2)<1e-9);
        C = find(vecnorm(nurbsStr(i).collocPts-[1,1], Inf, 2)<1e-9);
        D = find(vecnorm(nurbsStr(i).collocPts-[0,1], Inf, 2)<1e-9);
        if numel(B) > 0
            vecU(B*2-1:B*2)=PBC*[1; 0];
        end
        if numel(C) > 0
            vecU(C*2-1:C*2)=PBC*[1; 1];
        end
        if numel(D) > 0
            vecU(D*2-1:D*2)=PBC*[0; 1];
        end
    end
end

preU = unique(preU,'sorted');
preT = unique(preT,'sorted');
for i = 1:2
    U32{i} = unique(U32{i},'stable');
    T32{i} = unique(T32{i},'stable');
    U14{i} = unique(U14{i},'stable');
    T14{i} = unique(T14{i},'stable');

    temp = zeros(dem*length(U32{i}), 1);
    temp(1:2:end) = U32{i}*dem-1;
    temp(2:2:end) = U32{i}*dem;
    U32{i} = temp;

    temp = zeros(dem*length(T32{i}), 1);
    temp(1:2:end) = T32{i}*dem-1;
    temp(2:2:end) = T32{i}*dem;
    T32{i} = temp;

    temp = zeros(dem*length(U14{i}), 1);
    temp(1:2:end) = U14{i}*dem-1;
    temp(2:2:end) = U14{i}*dem;
    U14{i} = temp;

    temp = zeros(dem*length(T14{i}), 1);
    temp(1:2:end) = T14{i}*dem-1;
    temp(2:2:end) = T14{i}*dem;
    T14{i} = temp;
end
unU = setdiff((1:dem*numU)',preU);
unT = setdiff((1:dem*numT)',preT);
end