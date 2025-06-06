function [xi_u, flag] = JudgeSingular2D(nurbsStr,xiParam, range, ip, id, el)
% Function: JudgeSingular2D
% Description: Determines if a collocation point is at a singular location
%               within a 2D element.
%
% Input:
%   xiParam: Parameter coordinate of the collocation point.
%   range: Range of parameter coordinates for the current element.
%   ip: Index of the current boundary structure.
%   id: Index of the current element.
%   el: Index of element.
%
% Output:
%   xi_u: Parameter coordinate of the collocation point.
%   flag: Flag indicating singularity (0 for not singular, 1 for singular).

% Initialize
flag = 0;
xi_u = 0;  % A placeholder for regular integral.

% Check if the collocation point is within the parameter range of the element
if (xiParam(1) <= range(2)) && (xiParam(1) >= range(1) || ...
        (el == nurbsStr(id).numElem && xiParam(1) == 0)) && ip == id
    xi_u = xiParam(1);
    if el == nurbsStr(id).numElem && xiParam(1) == 0
        xi_u = 1;
    end
    flag = 1;
end

end
