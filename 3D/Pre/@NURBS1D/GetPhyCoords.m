
function PhyCoords = GetPhyCoords(obj,xiP)
%�������������ꡡ�����Ӧ����ռ������
%  xiP Ϊ����Ĳ�������
%  PhyCoords Ϊ����ռ������ֵ Ϊ������

% �ж������ĵ�Ԫ ���õ�Ԫ�Ĳ������귶Χ
L1 = logical(obj.elementVertex(:,1)-xiP<=0);
L2 = logical(obj.elementVertex(:,2)-xiP>=0);

e = find((L1 & L2),1);

gxmin = obj.elementVertex(e,1);
gxmax = obj.elementVertex(e,2);

% �ҳ���Ӧ�õ�Ԫ�ĸ�˹����
uu_hat = (2*xiP - gxmin - gxmax)/(gxmax-gxmin);

% �����Ӧ���κ��� ���õ�Ԫ�����Ŀ��Ƶ�
R = obj.Shapes(uu_hat,e);
cpts = obj.controlPts( obj.elementNode(e,:),:);

PhyCoords = R * cpts;

end
