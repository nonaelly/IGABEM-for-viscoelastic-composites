
function [R, dRdxi, dRdeta] = ShapesAnd1stDers(obj,xiG,etaG,ele)

% 给定单元号 和 高斯点坐标 返回非零形函数对参数坐标一阶导数的值
% xiG,etaG 为高斯点坐标
% ele为单元编号
% 输出 R dRdxi 和 dRdeta 为 行向量

xmin = obj.elementVertex(ele,1);
xmax = obj.elementVertex(ele,3);
ymin = obj.elementVertex(ele,2);
ymax = obj.elementVertex(ele,4);
    
indexU = find(ismember(obj.uKnot,xmin),1,'last') - 1;
indexV = find(ismember(obj.vKnot,ymin),1,'last') - 1;

Xi = 0.5 * ((xmax - xmin ) * xiG + xmax + xmin);
Eta = 0.5 * ((ymax - ymin ) * etaG + ymax + ymin);

[R, dRdxi, dRdeta] = NURBS2DBasisDers([Xi; Eta],obj.p,obj.q,obj.uKnot,obj.vKnot,obj.weights',indexU,indexV);

end