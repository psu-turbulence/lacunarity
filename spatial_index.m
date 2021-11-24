function [r_a, h] = spatial_index(R, L, tol)
% Cutoff point and heterogeneity index assignment
%
% INPUTS:
% R 			- Box values from lacunarity.m
% L             - Lacunarity values from lacunarity.m
% tol           - Tolerance for L = 1
%
% OUTPUTS:
% r_a 			- Cutoff point
% H             - Heterogeneity index [0, 1]
%
% CREATED:
% Ryan Scott
% 03/23/2021

d_L = gradient(L, R);
idx = squeeze(find(d_L > 0 & L/L(1) < 1 + tol, 1) - 1);
L = L/L(1);
if isempty(idx)
    r_a = R(end);
    l_w = (1/length(R))*sum(L.*R);
elseif (idx < 1)
    r_a = R(end);
    l_w = (1/length(R))*sum(L.*R);
else
    r_a = R(idx);
    l_w = (1/idx)*sum(L(1:idx).*R(1:idx));
end
h = 1 - (2*l_w)/(1 + r_a);

end
