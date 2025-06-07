
function PhyCoords = GetPhyCoords(obj,xiP,etaP)
%　给定参数坐标　求出对应物理空间的坐标
%  xiP etaP 为输入的参数坐标
%  PhyCoords 为物理空间的坐标值 为行向量

% 判断所属的单元 及该单元的参数坐标范围
L1 = logical(obj.elementVertex(:,1)-xiP<=0);
L2 = logical(obj.elementVertex(:,3)-xiP>=0);
L3 = logical(obj.elementVertex(:,2)-etaP<=0);
L4 = logical(obj.elementVertex(:,4)-etaP>=0);

e = find((L1 & L2 & L3 & L4),1);

gxmin = obj.elementVertex(e,1);
gymin = obj.elementVertex(e,2);

indexU = find(ismember(obj.uKnot,gxmin),1,'last') - 1;
indexV = find(ismember(obj.vKnot,gymin),1,'last') - 1;

% 计算对应的形函数 及该单元所属的控制点
R = NURBS2Dbasis([xiP; etaP],obj.p,obj.q,obj.uKnot,obj.vKnot,obj.weights',indexU,indexV);
cpts = obj.controlPts( obj.elementNode(e,:),:);

PhyCoords = R * cpts;

end
