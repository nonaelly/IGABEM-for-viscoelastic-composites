function DxDxiG = GetPhyToGaussDers(obj,xiG,ele)

% 输入 单元 编号 和 高斯点坐标
% 得到 物理空间坐标对 高斯空间坐标的 一阶导数 二阶导数
% 得到的结果 均为列向量 3 X 1

dRdxi = obj.Shape1stDers(xiG,ele);

xmin = obj.elementVertex(ele,1);
xmax = obj.elementVertex(ele,2);

cpts = obj.controlPts( obj.elementNode(ele,:),:);
% calculate the Jacobian of the transformation
dxdxi = dRdxi*cpts;  % 1 X 1 matrix
% dxdxi = [  dxdxi   dydxi  dzdxi

DxDxiG = dxdxi*(xmax-xmin)/2;

end
