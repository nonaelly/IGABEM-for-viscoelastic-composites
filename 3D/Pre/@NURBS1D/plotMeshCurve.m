function plotMeshCurve(obj,s1)

resolution = 90;
% discretize the xi and eta directions
xiVec   = linspace(0,max(obj.uKnot),resolution);

x1  = zeros(resolution,1);
y1  = zeros(resolution,1);
z1  = zeros(resolution,1);

% NURBS curves of knot lines corresponding to xi direction

projcoord = nurb2proj(obj.controlPts, obj.weights);
dim=size(projcoord,2);

for i=1:resolution
    
    xi = xiVec(i);
    tem  = CurvePoint(obj.noPtsU-1,obj.p,obj.uKnot,projcoord,dim,xi);
    x1(i) = tem(1)/tem(4);
    y1(i) = tem(2)/tem(4);
    z1(i) = tem(3)/tem(4);
end

plot3(x1,y1,z1,'r-','LineWidth',2);
hold on

if strcmp(s1,'PlotCPs')
    
    plot3(obj.controlPts(:,1),obj.controlPts(:,2),obj.controlPts(:,3),'bo',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',9,'LineWidth',2);
    
    hold on
    
    plot3(obj.controlPts(:,1),obj.controlPts(:,2),obj.controlPts(:,3),'b--','LineWidth',1);

end

% hold off
% axis square
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view(3)

end

function projcoord = nurb2proj( controlPoints, weights)
%--------------------------------------------------------------
% function projcoord = nurb2proj(nob, controlPoints, weights)
% transform NURBS data into projective coordinates
% INPUT:
% controlPoints: vector of control points (1 per row)
% weights :    : column vector of weights
% OUTPUT:
% projcoord    : matrix with projective coordinates
%--------------------------------------------------------------
projcoord = controlPoints;
for i=1:size(controlPoints,1)
    projcoord(i,:) = projcoord(i,:)*weights(i);
end
projcoord = [projcoord, weights];
end

function S = CurvePoint(n,p,U,P,dim,u)
%--------------------------------------------------------------
%function S = SolidPoint(n,p,U,m,q,V,P,u,v)
% this can be used for B-Spline solid
% or NURBS solid in projective coordinates (correct variable dim!)
%INPUT:
% n         : number ob basis functions -1 !  - x-direction
%        NURBS-Book: n+1 # basis, np max index (startindex 0)
%        here        n   # basis and max index (startindex 1)
% p          : degree of the basis functions - x-direction
% U          : knotvector - x-direction
% m          : number ob basis functions -1 !  - y-direction
% q          : degree of the basis functions - y-direction
% V          : knotvector - y-direction
% P          : control points
% dim        : dimension of control points
% u          : xi-coordinate
% v          : eta-coordinate
% w          : zeta-coordinate
%OUTPUT:
% S          : coordinates of the point on the solid
%--------------------------------------------------------------

% find spans of (u,v,w)
uspan = FindSpan(n,p,u,U);
% compute non-zero B-spline basis functions
Nu    = BasisFun(uspan,u,p,U);
% compute point on solid using B-spline interpolation
uind = uspan -p;
S    = zeros(1,dim);
for i=0:p
    CP   = P(uind+i+1,:);
    S    = S + Nu(i+1) * CP;
end
end

function N = BasisFun(i,u,p,U)
%--------------------------------------------------------------
% function N = BasisFun(i,p,u,U)
% NURBS-Book (algorithm A2.2)
% evalute nonzero basis functions
% INPUT:
% i          : current knotspan
% u          : evaluation point
% p          : degree of the basis functions
% U          : knot vector (row vector)
% OUTPUT:
% N          : row vector (dim p+1)
%              values of the basis function N_(i-p) ... N_(i)
%              at the evaluation point u
%--------------------------------------------------------------
N=zeros(1,p+1);
N(1)=1;
left=zeros(1,p+1);
right=zeros(1,p+1);
for j=1:p
    left(j+1) = u-U(i+1-j+1);
    right(j+1) = U(i+j+1)-u;
    saved = 0;
    for r=0:j-1
        temp = N(r+1)/(right(r+2)+left(j-r+1));
        N(r+1) = saved + right(r+2)*temp;
        saved = left(j-r+1)*temp;
    end
    N(j+1) = saved;
end
end

function knotSpanIndex = FindSpan(n,p,u,U)
%--------------------------------------------------------------
% function knotSpanIndex = FindSpan(n,p,u,U)
% NURBS-Book (algorithm A2.1)
% find the knot span index for one variable u
% INPUT:
% n          : number of basis function -1
%        NURBS-Book: np+1 # basis, np max index (startindex 0)
%        here        np   # basis and max index (startindex 1)
% p          : degree of the basis functions
% u          : evaluation point
% U          : knot vector (row vector)
% OUTPUT:
% knotSpanIndex : index of knot span
%--------------------------------------------------------------
if (u == U(n+2))
    knotSpanIndex= n;
    return
end
low = p;
high = n+1;
mid = floor((low + high)/2);
while (u <U(mid+1) || u >= U(mid+2) )
    if( u < U(mid+1))
        high = mid;
    else
        low = mid;
    end
    mid = floor((low+high)/2);
end
knotSpanIndex = mid;
end
