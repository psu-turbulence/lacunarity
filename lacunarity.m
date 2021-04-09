function [L, Ln, R, Rn, Z, Zn] = lacunarity(data, n_p)
% 1, 2, 3 dimensional lacunarity calculation
%
% INPUTS:
% data          - Input data, must be 2D matrix of element height
%
% Optional:
% n_p 			- Number of box sizes, must be lower than data size (default: n_p = max dimension of data)
%
% OUTPUTS:
% L 			- Vector of lacunarity values
% Ln            - Vector of lacunarity values normalized by L(1)
% R 			- Vector of box sizes
% Rn            - Vector of box sizes normalized by R(1)
% Z             - Matrix of first four moments
% Zn            - Normalized matrix of first four moments
%
% CREATED:
% Ryan Scott
% 03/23/2021
%
%EDITS:
% Sarah Smith 04/09/2021: 
%   - added 'nargin' check to account for full data rather than requre n_p box size designation
%%

s = size(data);
s_max = max(s) - 1;

if nargin<2
    n_p = s_max;
end
if (n_p > max(s))
    fprintf (2, 'Number of points must be less than data size.\n');
    return
end

R = unique(round(linspace(1, s_max, n_p)));
L = 1:length(R);
Z = NaN(length(R), 4);

if (length(s) == 3)
    
    % Flatten 3D data
    data = reshape(data, s(1), numel(data)/s(1));
    for b = 1:length(R)
        r = R(b);
        A = data;
        
        % Direction 1
        if (r < s(1))
            F = spdiags(repmat(ones(s(1), 1), 1, r), 0:(r - 1), s(1), s(1));
            A = F*A;
            A((end - r + 2):end, :) = [];
            A = reshape(A', s(2), numel(A)/s(2));
        else
            F = ones(s(1), s(1));
            A = F*A;
            A = reshape(A', s(2), numel(A)/s(2));
        end
        
        % Direction 2
        if (r < s(2))
            F = spdiags(repmat(ones(s(2), 1), 1, r), 0:(r - 1), s(2), s(2));
            A = F*A;
            A((end - r + 2):end, :) = [];
            A = reshape(A', s(3), numel(A)/s(3));
        else
            F = ones(s(2), s(2));
            A = F*A;
            A = reshape(A', s(3), numel(A)/s(3));
        end
        
        % Direction 3
        if (r < s(3))
            F = spdiags(repmat(ones(s(3), 1), 1, r), 0:(r - 1), s(3), s(3));
            A = F*A;
            A((end - r + 2):end, :) = [];
        else
            F = ones(s(3), s(3));
            A = F*A;
        end
        
        % Compute statistics
        z1 = mean(A(:)/r, 'omitnan');
        z2 = var(A(:)/r, 'omitnan');
%         z3 = skewness(A(:)/r, 'omitnan');
%         z4 = kurtosis(A(:)/r, 'omitnan');
        z3 = skewness(A(:)/r);      %'omitnan' inappropriate For tested version of MATLAB 2020a
        z4 = kurtosis(A(:)/r);
        Z(b, :) = [z1, z2, z3, z4];
        L(b) = 1 + z2/(z1^2);

        % Progress
        disp(r/s_max);
    
    end
elseif (length(s) == 2 && ~any(s == 1))
    for b = 1:length(R)
        r = R(b);
        A = data;
        
        % Direction 1
        if (r < s(1))
            F = spdiags(repmat(ones(s(1), 1), 1, r), 0:(r - 1), s(1), s(1));
            A = F*A;
            A((end - r + 2):end, :) = [];
            A = A';
        else
            F = ones(s(1), s(1));
            A = F*A;
            A = A';
        end
        
        % Direction 2
        if (r < s(2))
            F = spdiags(repmat(ones(s(2), 1), 1, r), 0:(r - 1), s(2), s(2));
            A = F*A;
            A((end - r + 2):end, :) = [];
        else
            F = ones(s(2), s(2));
            A = F*A;
        end
        
        % Compute statistics
        z1 = mean(A(:)/r, 'omitnan');
        z2 = var(A(:)/r, 'omitnan');
        z3 = skewness(A(:)/r, 'omitnan');
        z4 = kurtosis(A(:)/r, 'omitnan');
        Z(b, :) = [z1, z2, z3, z4];
        L(b) = 1 + z2/(z1^2);

        % Progress
        disp(r/s_max);
    end
elseif (length(s) == 2 && any(s == 1))
    s_1D = max(s);
    if (s(1) == 1)
        data = data';
    end
    for b = 1:length(R)
        r = R(b);
        A = data;
        
        F = spdiags(repmat(ones(s_1D, 1), 1, r), 0:(r - 1), s_1D, s_1D);
        A = F*A;
        A((end - r + 2):end, :) = [];
        
        % Compute statistics
        z1 = mean(A(:)/r, 'omitnan');
        z2 = var(A(:)/r, 'omitnan');
%         z3 = skewness(A(:)/r, 'omitnan');
%         z4 = kurtosis(A(:)/r, 'omitnan');
        z3 = skewness(A(:)/r);
        z4 = kurtosis(A(:)/r);
        Z(b, :) = [z1, z2, z3, z4];
        L(b) = 1 + z2/(z1^2);

        % Progress
        disp(r/s_max);
    end
else
    fprintf (2, 'Requires 1D 2D or 3D data.\n');
    return
end

% Normalize
Ln = L/L(1);
Rn = (R - R(1))/(s_max - R(1));
Zn = [Z(:, 1)/Z(1, 1),Z(:, 2)/Z(1, 2),Z(:, 3)/Z(1, 3),Z(:, 4)/Z(1, 4)];

end