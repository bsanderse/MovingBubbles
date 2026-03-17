function phi = wendland(r, k)
% WENDLAND radial basis functions
% r : scaled distance (||x-x_i|| / epsilon)
% k : smoothness parameter (0,1,2)
%
% Output:
% phi : Wendland basis value

r = max(r,0);
phi = zeros(size(r));

switch k
    
    % C^0
    case 0
        mask = r < 1;
        phi(mask) = (1 - r(mask)).^2;
        
    % C^2
    case 1
        mask = r < 1;
        phi(mask) = (1 - r(mask)).^4 .* (4*r(mask) + 1);
        
    % C^4
    case 2
        mask = r < 1;
        phi(mask) = (1 - r(mask)).^6 .* (35*r(mask).^2 + 18*r(mask) + 3) / 3;
        
    otherwise
        error('k must be 0,1,2')
end
end