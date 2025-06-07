function [CollocPattern, sizeBasis] = GluePatches(patchBoundaries,shellNurbs)

% patchBoundaries 的格式为
% patchAList patchB edgeAList edgeBList
% 将 patchB 与 patchAList 粘接起来
% 一个片一个片的粘
% 要求粘接处 参数方向一致 ##########################

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
    
    % 获取片和边的信息
    patchAList = patchBoundaries{boundaryIndex,1}; 
    patchB     = patchBoundaries{boundaryIndex,2};  
    edgeAList  = patchBoundaries{boundaryIndex,3}; 
    edgeBList  = patchBoundaries{boundaryIndex,4}; 
    DirList    = patchBoundaries{boundaryIndex,5};
    
    nodesPattern = zeros(1, shellNurbs(patchB).noPts); 
    % 改变 patchB 的控制点编号

    nodesA = [];
    nodesB = [];
    patchesSeen = union(patchesSeen, patchB); % 取并集
    
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
    
    % 将片与片之间连接处的nodesA和nodesB按相同顺序(基于patchB的顺序)排序
    [nodesB,sI] = sort(nodesB);
    nodesA = nodesA(sI);
    
    curBdryNode = 0;
    for nodeIndex=1:length(nodesA)
        % shift the basis functions indices in nodesPattern
        % nodesB里与nodesA里对应位置的自由度与A相同，其余的为当前自由度加上dimbasisA的
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
    
    % 将A patch 与B patch 控制点编号统一起来 ,注意两片的边界上的控制点是共享的,并且以A patch的控制点为准
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

