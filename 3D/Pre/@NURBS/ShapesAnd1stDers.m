
function [R, dRdxi, dRdeta] = ShapesAnd1stDers(obj,xiG,etaG,ele)

% ������Ԫ�� �� ��˹������ ���ط����κ����Բ�������һ�׵�����ֵ
% xiG,etaG Ϊ��˹������
% eleΪ��Ԫ���
% ��� R dRdxi �� dRdeta Ϊ ������

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