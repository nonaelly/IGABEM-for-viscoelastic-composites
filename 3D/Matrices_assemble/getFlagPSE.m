function flag = getFlagPSE(xi_u, xi_v, range)

flag = [0,0,0,0];
tol = 1e-10;
index = [4,2,1,3];
for i = 1:2
    if abs(xi_u - range(i)) < tol
        flag(index(i)) = 1;
    end
    if abs(xi_v - range(i+2)) < tol
        flag(index(i+2)) = 1;
    end
end


