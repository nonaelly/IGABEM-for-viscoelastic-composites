
function [R, dRdxi, dRdeta, dR2dxi2, dR2deta2, dR2dxieta] = Shapes1st2ndDers(obj,xiG,etaG,ele)

% ������Ԫ�� �� ��˹������ ���ط����κ���,�Բ�������һ�׵���,�Բ����Ķ��׵�����ֵ
% xiG,etaG Ϊ��˹������
% eleΪ��Ԫ���
% ��� �� ��Ϊ ������

xmin = obj.elementVertex(ele,1);
xmax = obj.elementVertex(ele,3);
ymin = obj.elementVertex(ele,2);
ymax = obj.elementVertex(ele,4);
    
indexU = find(ismember(obj.uKnot,xmin),1,'last') - 1;
indexV = find(ismember(obj.vKnot,ymin),1,'last') - 1;

Xi = 0.5 * ((xmax - xmin ) * xiG + xmax + xmin);
Eta = 0.5 * ((ymax - ymin ) * etaG + ymax + ymin);

[R, dRdxi, dRdeta, dR2dxi2, dR2deta2, dR2dxieta] = NURBS2DBasis2ndDers([Xi; Eta],obj.p,obj.q,obj.uKnot,obj.vKnot,obj.weights',indexU,indexV);

% notify(obj,'ShapeUsedAgain');

end
