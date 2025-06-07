
function PhyCoords = GetPhyCoords(obj,xiP)
%　给定参数坐标　求出对应物理空间的坐标
%  xiP 为输入的参数坐标
%  PhyCoords 为物理空间的坐标值 为行向量

% 判断所属的单元 及该单元的参数坐标范围
L1 = logical(obj.elementVertex(:,1)-xiP<=0);
L2 = logical(obj.elementVertex(:,2)-xiP>=0);

e = find((L1 & L2),1);

gxmin = obj.elementVertex(e,1);
gxmax = obj.elementVertex(e,2);

% 找出对应该单元的高斯坐标
uu_hat = (2*xiP - gxmin - gxmax)/(gxmax-gxmin);

% 计算对应的形函数 及该单元所属的控制点
R = obj.Shapes(uu_hat,e);
cpts = obj.controlPts( obj.elementNode(e,:),:);

PhyCoords = R * cpts;

end
