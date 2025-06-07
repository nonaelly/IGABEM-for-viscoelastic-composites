function [CollocPattern, sizeBasis] = GluePatches(patchBoundaries,shellNurbs)

% patchBoundaries �ĸ�ʽΪ
% patchAList patchB edgeAList edgeBList
% �� patchB �� patchAList ճ������
% һ��Ƭһ��Ƭ��ճ
% Ҫ��ճ�Ӵ� ��������һ�� ##########################

% patchBoundaries = {1,2,2,4;
%     2,3,2,2;
%     [1,3],4,[4,4],[4,2];
%     [1,2,3,4],5,[1,1,1,1],[4,3,2,1];
%     [1,2,3,4],6,[3,3,3,3],[4,3,2,1]};

numBoundaries = size(patchBoundaries,1);

curShift = shellNurbs(1).noPts;
overlapCounter = 0;
patchesSeen = [];

CollocPattern = cell(length(shellNurbs),1);
CollocPattern{1} = 1:shellNurbs(1).noPts;
for boundaryIndex = 1:numBoundaries 
    
    % ��ȡƬ�ͱߵ���Ϣ
    patchAList = patchBoundaries{boundaryIndex,1}; 
    patchB     = patchBoundaries{boundaryIndex,2};  
    edgeAList  = patchBoundaries{boundaryIndex,3}; 
    edgeBList  = patchBoundaries{boundaryIndex,4}; 
    DirList    = patchBoundaries{boundaryIndex,5};
    
    nodesPattern = zeros(1, shellNurbs(patchB).noPts); 
    % �ı� patchB �Ŀ��Ƶ���

    nodesA = [];
    nodesB = [];
    patchesSeen = union(patchesSeen, patchB); % ȡ����
    
    for indexPatch = 1:length(patchAList) 
        
        patchA = patchAList(indexPatch);
        edgeA  =  edgeAList(indexPatch); 
        edgeB  =  edgeBList(indexPatch); 
        FlagDir =    DirList(indexPatch); 
        
        patchesSeen = union(patchesSeen, patchA);
        
        [nodesAadd, ~] = shellNurbs(patchA).sortNodesAndElemsOfEdge(edgeA,1);
        [nodesBadd, ~] = shellNurbs(patchB).sortNodesAndElemsOfEdge(edgeB,FlagDir);
        
        nodesA = [nodesA; nodesAadd];
        nodesB = [nodesB; nodesBadd];

    end
    
    % ��Ƭ��Ƭ֮�����Ӵ���nodesA��nodesB����ͬ˳��(����patchB��˳��)����
    [nodesB,sI] = sort(nodesB);
    nodesA = nodesA(sI);
    
    curBdryNode = 0;
    for nodeIndex=1:length(nodesA)
        % shift the basis functions indices in nodesPattern
        % nodesB����nodesA���Ӧλ�õ����ɶ���A��ͬ�������Ϊ��ǰ���ɶȼ���dimbasisA��
        prevBdryNode = curBdryNode;
        curBdryNode = nodesB(nodeIndex);
        nodesPattern(prevBdryNode+1:curBdryNode-1) = ((prevBdryNode+1):(curBdryNode-1)) + curShift;
        nodesPattern(curBdryNode) = nodesA(nodeIndex);
        if prevBdryNode < curBdryNode
            curShift = curShift - 1;
        end
    end
    % shift the indices after the last boundary node
    nodesPattern(curBdryNode+1:end) = ((curBdryNode+1):shellNurbs(patchB).noPts) + curShift;
    
    % ��A patch ��B patch ���Ƶ���ͳһ���� ,ע����Ƭ�ı߽��ϵĿ��Ƶ��ǹ����,������A patch�Ŀ��Ƶ�Ϊ׼
    % update the nodesGlobal in patchB according to nodesPattern    
    
    shellNurbs(patchB).eleNodeGlobal = nodesPattern(shellNurbs(patchB).elementNode);
    CollocPattern{patchB} = nodesPattern;
    overlapCounter = overlapCounter + length(unique(nodesA));
    
    dimBasis = 0;
    for iobj = 1:length(patchesSeen)
        dimBasis = dimBasis + shellNurbs(iobj).noPts;
    end
    curShift = dimBasis - overlapCounter;
    
end

dimBasis = 0;
for iobj = 1:length(shellNurbs)
    dimBasis = dimBasis + shellNurbs(iobj).noPts;
end
sizeBasis = dimBasis - overlapCounter;

end

