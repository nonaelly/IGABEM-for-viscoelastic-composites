
function PhyCoords = GetPhyCoords(obj,xiP,etaP)
%�������������ꡡ�����Ӧ����ռ������
%  xiP etaP Ϊ����Ĳ�������
%  PhyCoords Ϊ����ռ������ֵ Ϊ������

% �ж������ĵ�Ԫ ���õ�Ԫ�Ĳ������귶Χ
L1 = logical(obj.elementVertex(:,1)-xiP<=0);
L2 = logical(obj.elementVertex(:,3)-xiP>=0);
L3 = logical(obj.elementVertex(:,2)-etaP<=0);
L4 = logical(obj.elementVertex(:,4)-etaP>=0);

e = find((L1 & L2 & L3 & L4),1);

gxmin = obj.elementVertex(e,1);
gymin = obj.elementVertex(e,2);

indexU = find(ismember(obj.uKnot,gxmin),1,'last') - 1;
indexV = find(ismember(obj.vKnot,gymin),1,'last') - 1;

% �����Ӧ���κ��� ���õ�Ԫ�����Ŀ��Ƶ�
R = NURBS2Dbasis([xiP; etaP],obj.p,obj.q,obj.uKnot,obj.vKnot,obj.weights',indexU,indexV);
cpts = obj.controlPts( obj.elementNode(e,:),:);

PhyCoords = R * cpts;

end
