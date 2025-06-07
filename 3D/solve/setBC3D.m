function [preU, preT, unU, unT, vecT, vecU, tranU, tranT, couUM, couTM, couT1M...
    , couUS, couTS, couT1S] = setBC3D(BC,nurbsStr,pressure,dirichlet,collocNew)
% Generate boundary condition for each 2D boundary.
D = 3;
[indU, indT, indNumU, indNumT, indNumB] = initializeIndices(nurbsStr, collocNew);
numU = indNumU(end);
numT = indNumT(end);
numBou = length(collocNew);
vecU = zeros(D*numU,1);
vecT = zeros(D*numT,1);
preU = [];
preT = [];
indDir = 0;
indPre = 0;
tranU = eye(D*numU);
tranT = eye(D*numT);
% patchLink : v=0,u=1,v=1,u=0
[couUM, couTM, couT1M, couUS, couTS, couT1S]=deal([]);
for i = 1 : length(nurbsStr)

    uConn = nurbsStr(i).globU(unique(reshape(nurbsStr(i).uConn, 1, [])))';
    tConn = nurbsStr(i).globT(unique(reshape(nurbsStr(i).tConn, 1, [])))';

    switch BC(i)
        case {'x','dx'}
            tempU = uConn*D - D + 1;
            preU = [preU; tempU];
            tempT = [tConn*D; tConn*D - D + 2];
            preT = [preT; tempT];
            if strcmpi(BC(i),'x')
                indDir = indDir + 1;
                vecU(tempU) = dirichlet(indDir);
            end
        case {'y','dy'}
            tempU = uConn*D - D + 2;
            preU = [preU; tempU];
            tempT = [tConn*D; tConn*D - D + 1];
            preT = [preT; tempT];
            if strcmpi(BC(i),'y')
                indDir = indDir + 1;
                vecU(tempU) = dirichlet(indDir);
            end
        case {'z','dz'}
            tempU = uConn*D;
            preU = [preU; tempU];
            tempT = [tConn*D - D + 2; tConn*D - D + 1];
            preT = [preT; tempT];
            if strcmpi(BC(i),'z')
                indDir = indDir + 1;
                vecU(tempU) = dirichlet(indDir);
            end
        case {'f','px','py','pz','pr'}
            tempT = [tConn*D-2; tConn*D-1; tConn*D];
            preT = [preT; tempT];
            switch BC(i)
                case 'px'
                    indPre = indPre + 1;
                    tempP = tConn;
                    vecT(tempP*D-2) = -pressure(indPre);
                case 'py'
                    indPre = indPre + 1;
                    tempP = tConn;
                    vecT(tempP*D-1) = -pressure(indPre);
                case 'pz'
                    indPre = indPre + 1;
                    tempP = tConn;
                    vecT(tempP*D) = -pressure(indPre);
            end
            if strcmpi(BC(i),'pr')
                indPre = indPre + 1;
                tempP = tConn;
                normP = nurbsStr(i).normCollocT(unique(reshape(nurbsStr(i).tConn, 1, [])),:);
                vecT(tempP*D-2) = -pressure(indPre) * normP(:,1);
                vecT(tempP*D-1) = -pressure(indPre) * normP(:,2);
                vecT(tempP*D) = -pressure(indPre) * normP(:,3);
            end
        case 'd' % Fix
            tempU = [uConn*D-2; uConn*D-1; uConn*D];
            [preU] = [preU; reshape(tempU,[],1)];
        case 'dr' % Skew boundary
            % Similar as 'dx'
            for j = 1 : length(uConn)
                cosT = nurbsStr(i).normCollocT(j, 1);
                sinT = nurbsStr(i).normCollocT(j, 2);
                R_inv = [cosT -sinT 0 ;sinT cosT 0; 0 0 1];
                idx = [uConn(j)*D-2; uConn(j)*D-1; uConn(j)*D];
                tranU(idx, idx) = R_inv;

                idx = [tConn(j)*D-2; tConn(j)*D-1; tConn(j)*D];
                tranT(idx, idx) = R_inv;
            end
            tempU = uConn*D-2;
            [preU] = [preU; reshape(tempU,[],1)];
            tempT = [tConn*D; tConn*D - D + 2];
            [preT] = [preT; tempT];
        case {'cm', 'cs'} % couple faces (master & slave)
            switch BC(i)
                case 'cm'
                    % uS=uM, tS=-tM
                    tempU = [uConn*D-2; uConn*D-1; uConn*D];
                    couUM = reshape(tempU,[],1);
                    tempT = [tConn*D-2; tConn*D-1; tConn*D];
                    couTM = tempT;

                    % We should also remember the lines where uS=uM, tS=tM
                    couT1M = [];
                    for j=1:4
                        k = nurbsStr(i).patchLink(j);
                        nU = (size(nurbsStr(k).knotU, 2)-nurbsStr(k).p-1);
                        nV = (size(nurbsStr(k).knotV, 2)-nurbsStr(k).q-1);

                        % find the place of the i-th patch in k-th patch.
                        switch find(nurbsStr(k).patchLink==i)
                            case 1 % v=0
                                tCTemp = nurbsStr(k).t2uConn(1:nU);
                            case 2 % u=1
                                tCTemp = nurbsStr(k).t2uConn(nU:nU:end);
                            case 3 % v=1
                                tCTemp = nurbsStr(k).t2uConn(end+1-nU:end);
                            case 4 % u=0
                                tCTemp = nurbsStr(k).t2uConn(1:nU:end);
                        end
                        tCTemp = nurbsStr(k).globT(tCTemp)';
                        tempT = [tCTemp*D-2; tCTemp*D-1; tCTemp*D];
                        [couT1M] = [couT1M; tempT];
                    end

                case 'cs'
                    uConnFlip = [];
                    tConnFlip = [];
                    % we should flip based on the face (1-4 / 5,6)
                    % uS=uM, tS=-tM
                    if ismember(i,1:4)
                        nU = (size(nurbsStr(i).knotU, 2)-nurbsStr(i).p-1);
                        nV = (size(nurbsStr(i).knotV, 2)-nurbsStr(i).q-1);
                        for j=1:nV
                            [uConnFlip] = [uConnFlip; flip(uConn(1+nU*(j-1) : nU*j))];
                            [tConnFlip] = [tConnFlip; flip(tConn(1+nU*(j-1) : nU*j))];
                        end
                        tempU = [uConnFlip*D-2; uConnFlip*D-1; uConnFlip*D];
                        couUS = reshape(tempU,[],1);
                        tempT = [tConnFlip*D-2; tConnFlip*D-1; tConnFlip*D];
                        couTS = tempT;

                        % We should also remember the lines where uS=uM, tS=tM
                        for j=[1,4,3,2]
                            k = nurbsStr(i).patchLink(j);
                            nU = (size(nurbsStr(k).knotU, 2)-nurbsStr(k).p-1);
                            nV = (size(nurbsStr(k).knotV, 2)-nurbsStr(k).q-1);

                            % find the place of the i-th patch in k-th patch.
                            switch find(nurbsStr(k).patchLink==i)
                                case 1 % v=0
                                    tCTemp = nurbsStr(k).t2uConn(1:nU);
                                case 2 % u=1
                                    tCTemp = nurbsStr(k).t2uConn(nU:nU:end);
                                case 3 % v=1
                                    tCTemp = nurbsStr(k).t2uConn(end+1-nU:end);
                                case 4 % u=0
                                    tCTemp = nurbsStr(k).t2uConn(1:nU:end);
                            end
                            tCTemp = nurbsStr(k).globT(tCTemp)';
                            tempT = [tCTemp*D-2; tCTemp*D-1; tCTemp*D];
                            [couT1S] = [couT1S; tempT];
                        end

                    elseif ismember(i,5:6)
                        nU = (size(nurbsStr(i).knotU, 2)-nurbsStr(i).p-1);
                        nV = (size(nurbsStr(i).knotV, 2)-nurbsStr(i).q-1);
                        for j=1:nV
                            [uConnFlip] = [uConnFlip; flip(uConn(1+nU*(nV-j) : nU*(nV-j+1)))];
                            [tConnFlip] = [tConnFlip; flip(tConn(1+nU*(nV-j) : nU*(nV-j+1)))];
                        end
                        tempU = [uConnFlip*D-2; uConnFlip*D-1; uConnFlip*D];
                        couUS = reshape(tempU,[],1);
                        tempT = [tConnFlip*D-2; tConnFlip*D-1; tConnFlip*D];
                        couTS = tempT;

                        % We should also remember the lines where uS=uM, tS=tM
                        for j=[3,2,1,4]
                            k = nurbsStr(i).patchLink(j);
                            nU = (size(nurbsStr(k).knotU, 2)-nurbsStr(k).p-1);
                            nV = (size(nurbsStr(k).knotV, 2)-nurbsStr(k).q-1);

                            % find the place of the i-th patch in k-th patch.
                            switch find(nurbsStr(k).patchLink==i)
                                case 1 % v=0
                                    tCTemp = nurbsStr(k).t2uConn(1:nU);
                                case 2 % u=1
                                    tCTemp = nurbsStr(k).t2uConn(nU:nU:end);
                                case 3 % v=1
                                    tCTemp = nurbsStr(k).t2uConn(end+1-nU:end);
                                case 4 % u=0
                                    tCTemp = nurbsStr(k).t2uConn(1:nU:end);
                            end
                            tCTemp = nurbsStr(k).globT(tCTemp)';
                            tempT = [tCTemp*D-2; tCTemp*D-1; tCTemp*D];
                            [couT1S] = [couT1S; tempT];
                        end

                    else
                        error('over 6 patch')
                    end
            end
        otherwise
            error('Wrong BC')
    end
end


for i = 1 : length(nurbsStr)
    switch BC(i)% couple faces (master & slave)
        case 'cm'
        case 'cs'
            % add the preU on the poblic line
            preU = union(preU, couUS);
            preT = union(preT, couTS);
            preT = union(preT, couT1S);
    end
end

preU = unique(preU,'sorted');
preT = unique(preT,'sorted');
unU = setdiff((1:D*numU)',preU);
unT = setdiff((1:D*numT)',preT);
end