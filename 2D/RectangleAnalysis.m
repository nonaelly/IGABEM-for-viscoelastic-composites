% =========================================================================
% Analytical solution for a rectangular plate under uniaxial tension.
% The solution is based on a viscoelastic shear relaxation model:
%   G(t)=G_inf+G_1*exp(-t/tao)
%
% Only the displacement in the x-direction is considered.
%
% Input:
%   point : [x,y] coordinate to evaluate displacement
%   a     : length of plate (not directly used)
%   b     : height of plate (not directly used)
%   p     : applied pressure
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
function disp=RectangleAnalysis(point,a,b,p,GInf,G1,nu,t,tao)

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
    term1=(1/(2*q0G))*(1-exp(-q0G/q1G*ti));
    term2=(p1G/(2*q1G))*exp(-q0G/q1G*ti);
    term3=(1/(6*K+q0G))*(1-exp(-(6*K+q0G)/(6*K*p1G+q1G)*ti));
    term4=(p1G/(6*K*p1G+q1G))*exp(-(6*K+q0G)/(6*K*p1G+q1G)*ti);

    f(i)=term1+term2+3/2*(term3+term4);

    disp(i,1)=-f(i)*x*p;
    disp(i,2)=0;   % y-displacement is ignored
end

end
