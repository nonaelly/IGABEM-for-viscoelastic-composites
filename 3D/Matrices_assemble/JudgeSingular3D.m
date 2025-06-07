function [xi_u, xi_v, flag] = JudgeSingular3D(xiParam, range, ip, id, bouRange, flagLine, xiLine)

% On the same boundary.
flag1 = (id>=bouRange(1) && id<=bouRange(2));
% Two cases.
% 1: On the same patch
flag2 = (range(2) >= xiParam(1) && xiParam(1) >=range(1) && ...
    (range(4) >= xiParam(2) && xiParam(2) >= range(3)) && ip == id);
% 2: Not on the same patch, but in a common line
flag3 = flagLine && (range(2) >= xiLine(1) && xiLine(1) >=range(1) && ...
    (range(4) >= xiLine(2) && xiLine(2) >= range(3))) && id ~= ip;
flag = flag1 && (flag2||flag3);
if flag3
    xi_u = xiLine(1);
    xi_v = xiLine(2);
else
    xi_u = xiParam(1);
    xi_v = xiParam(2);
end

end