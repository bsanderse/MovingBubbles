function [r,drdt] = ellipse(theta,theta_dot,a,b)

    
    if (nargin<3)
        a = 2;
    end
    if (nargin<4)
        b = 1;
    end

    if (~exist('theta_dot','var'))
        theta_dot = 0;
    end


    cos_term = cos(theta);
    sin_term = sin(theta);
    
    % radius at each theta
    r = 1 ./ sqrt((cos_term.^2)/a^2 + (sin_term.^2)/b^2);
    
    % dr/dt using formula: dr/dt = -0.5 * r^3 * phi_dot * sin(2*(theta-phi0)) * (1/b^2 - 1/a^2)
    drdt = 0.5 * r.^3 * theta_dot .* sin(2*(theta)) * (1/b^2 - 1/a^2);

end
