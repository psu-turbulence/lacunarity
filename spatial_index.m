function [r_a, h] = spatial_index(L, R, tol)
% Cutoff point and heterogeneity index assignment
%
% INPUTS:
% L             - Lacunarity values from lacunarity.m
% R 			- Box values from lacunarity.m
% tol           - Gradient tolerance for finding dL/dR=0
%
% OUTPUTS:
% r_a 			- Cutoff point
% H             - Heterogeneity index [0, 1]
%
% CREATED:
% Ryan Scott
% 03/23/2021

d_L = abs(gradient(L3D3, R3D3));
idx = find(d_L < tol, 1);
if (isempty(idx))
    r_a = R(end);
    l_h = (1/length(R))*sum(L.*R);
else
    r_a = R3D3(idx);
    l_h = (1/idx)*sum(L(1:idx).*R(1:idx));
end
h = 1 - (2*l_h)/(1 + r_a);

end

