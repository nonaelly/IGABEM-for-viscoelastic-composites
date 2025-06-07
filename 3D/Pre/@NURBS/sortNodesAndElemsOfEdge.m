function [nodes, elements] = sortNodesAndElemsOfEdge(obj,noEdge,FlagDir)
% ѡȡ������ ���ղ��������� �ĵ�Ԫ�� �� ���Ƶ���
switch noEdge
    case 1
        elements = find(obj.elementLocation(:,1)==0);
        nodes = unique(obj.eleNodeGlobal(elements,obj.down_nodes)','stable');
    case 2
        elements = find(obj.elementLocation(:,2)==0);
        nodes = unique(obj.eleNodeGlobal(elements,obj.right_nodes)','stable');
    case 3
        elements = find(obj.elementLocation(:,3)==0);
        nodes = unique(obj.eleNodeGlobal(elements,obj.up_nodes)','stable');
    case 4
        elements = find(obj.elementLocation(:,4)==0);
        nodes = unique(obj.eleNodeGlobal(elements,obj.left_nodes)','stable');
end
if FlagDir == -1
    nodes = flip(nodes);
end

end
