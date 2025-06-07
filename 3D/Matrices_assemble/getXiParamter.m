function [Xi_u, Xi_v] = getXiParamter(xi_u, xi_v, range)

tol = 1e-15;
Xi_u = xi_u;
Xi_v = xi_v;

if xi_u == 0
    Xi_u = tol;
end
if xi_u == range(1) 
    Xi_u = range(1) + tol;
end
if xi_u == range(2)
    Xi_u = range(2) - tol;
end
if xi_v == range(3)
    Xi_v = range(3) + tol;
end
if xi_v == range(4)
    Xi_v = range(4) - tol;
end

end