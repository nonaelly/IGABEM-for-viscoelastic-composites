
function [dRdxi, J1] = Shape1stDersAndDRdx(obj,xiG,ele)

% 给定单元号 和 高斯点坐标 返回非零形函数对参数坐标一阶导数的值
% xiG 为高斯点坐标
% ele为单元编号
% 输出 dRdxi dRdeta 和 dRdzeta 为 行向量

ncpt = (obj.p + 1) ;

B    = zeros(ncpt,1);
dB = zeros(ncpt,1);

for i=1:obj.p+1
    B(i)=bernstein(obj.p,i,xiG);
    dB(i)=0.5*obj.p*(bernstein(obj.p-1,i-1,xiG)-bernstein(obj.p-1,i,xiG));
end

% form first ders basis functions using the Bezier extraction operator
B = (obj.C(:,:,ele)*B)';              % 1 x 16  CeBe
dBxi = (obj.C(:,:,ele)*dB(:,1))';     % 1 x 16     CeBexi

dN = zeros(1, ncpt);

% Multiply each B-spline function with corresponding weight
wgts = obj.weights(obj.elementNode(ele,:));
N = B.* wgts';    %  WeCeBe
w_sum = sum(N);   % sum W(xi，eta)

dN(1,:) = dBxi .* wgts';                    %  WeCeBexi

dw_xi = sum(dN(1,:));                       % sum W(dxi，eta, zeta)

dRdxi = dN(1,:)/w_sum - N*dw_xi/w_sum^2;    % dRdxi

% the first derivatives of shape functions w.r.t. parameter coordinates
xmin = obj.elementVertex(ele,1);
xmax = obj.elementVertex(ele,2);

dRdxi = dRdxi*2/(xmax-xmin);

cpts = obj.controlPts( obj.elementNode(ele,:),:);
% calculate the Jacobian of the transformation
dxdxi = dRdxi*cpts;  % 1 X 3 matrix
% dxdxi = [  dxdxi   dydxi  dzdxi ];
J1 = norm(dxdxi);

end

function B=bernstein(p,a,xi)

if p==0 && a==1
    B=1;
elseif p==0 && a~=1
    B=0;
else
    if a<1 || a>p+1
        B=0;
    else
        B1=bernstein(p-1,a,xi);
        B2=bernstein(p-1,a-1,xi);
        B=0.5*(1-xi)*B1+0.5*(1+xi)*B2;
    end
end
end
