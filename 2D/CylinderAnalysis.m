% =========================================================================
% Analytical solution of a thick-walled cylinder under internal pressure.
% The solution is based on a viscoelastic shear relaxation model:
%   G(t)=G_inf+G_1*exp(-t/tao)
%
% Input:
%   point : [x,y] coordinate to evaluate displacement
%   a     : inner radius
%   b     : outer radius
%   p     : internal pressure
%   GInf  : long-term shear modulus
%   G1    : decaying shear modulus
%   nu    : Poisson's ratio
%   t     : time array
%   tao   : relaxation time
%
% Output:
%   disp  : displacement [u_x,u_y] at point over time
%
% Author: [Wang Zhetong]
% =========================================================================
function disp=CylinderAnalysis(point,a,b,p,GInf,G1,nu,t,tao)

x=point(1);
y=point(2);
r=sqrt(x^2+y^2);
K=(GInf+G1)*2*(1+nu)/(3*(1-2*nu));

% Shear relaxation parameters
p1G=tao;
q0G=2*GInf;
q1G=2*(GInf+G1)*tao;
p0G=1;

disp=zeros(length(t),2);
f=zeros(length(t),1);

for i=1:length(t)
    ti=t(i);
    term1=3*r/(6*K+q0G);
    term2=(p1G/(q1G+6*K*p1G)-1/(6*K+q0G))*exp(-(6*K+q0G)/(q1G+6*K*p1G)*ti);
    term3=b^2/(q0G*r);
    term4=(q0G/q1G*p1G-1)*exp(-q0G/q1G*ti);

    f(i)=a^2*p/(b^2-a^2)*(term1+term2+term3*(1+term4));

    disp(i,1)=-f(i)*x/r;
    disp(i,2)=-f(i)*y/r;
end

end
