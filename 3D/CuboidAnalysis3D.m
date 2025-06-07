function disp = CuboidAnalysis3D(point, L, p, GInf, G1, nu, t, tao)
x = point(1);
y = point(2);
z = point(3);
r = sqrt(x^2+y^2);
K = (GInf+G1)*2*(1+nu)/(3*(1-2*nu));
a = L(1);
b = L(2);
c = L(3);

% G = G_Inf + G1*exp(-t/tao)
p1G = tao;
q0G = 2*GInf;
q1G = 2*(GInf+G1)*tao;
p0G = 1;
disp = zeros(length(t),1);
f = zeros(length(t),1);
for i=1:length(t)
    f(i)=((1/(9*K) + 2*p1G/3/q1G) + 2/3*(q1G-q0G*p1G)/(q1G*q0G)*(1-exp(-q0G/q1G*t(i))));
    disp(i,1) = -f(i)*p*a;
end
end