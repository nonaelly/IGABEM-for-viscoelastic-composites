function DxDxiG = GetPhyToGaussDers(obj,xiG,ele)

% ���� ��Ԫ ��� �� ��˹������
% �õ� ����ռ������ ��˹�ռ������ һ�׵��� ���׵���
% �õ��Ľ�� ��Ϊ������ 3 X 1

dRdxi = obj.Shape1stDers(xiG,ele);

xmin = obj.elementVertex(ele,1);
xmax = obj.elementVertex(ele,2);

cpts = obj.controlPts( obj.elementNode(ele,:),:);
% calculate the Jacobian of the transformation
dxdxi = dRdxi*cpts;  % 1 X 1 matrix
% dxdxi = [  dxdxi   dydxi  dzdxi

DxDxiG = dxdxi*(xmax-xmin)/2;

end
