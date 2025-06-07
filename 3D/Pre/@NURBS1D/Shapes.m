
function R = Shapes(obj,xiG,ele)

% 给定单元号 和 高斯点坐标 返回非零形函数的值
% xiG 为高斯点坐标
% ele为单元编号
% 输出 R 为 行向量

ncpt = (obj.p + 1) ;
B    = zeros(ncpt,1);
for i=1:obj.p+1
    B(i)=bernstein(obj.p,i,xiG);
end

% form B-spline basis functions using the Bezier extraction operator
B = (obj.C(:,:,ele)*B)';  % 1 x 16  CeBe

% Multiply each B-spline function with corresponding weight
wgts = obj.weights(obj.elementNode(ele,:));
R = B.* wgts';    %  WeCeBe
w_sum = sum(R);   % sum W(xi，eta)
R = R/w_sum;      % R = WeCeBe/sum W(xi，eta)

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


