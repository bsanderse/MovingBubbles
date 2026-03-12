function [x_new,dx_dt] = transform(x,t,mode)
%% choose transformation
% mode = 'rotation';
% options:
% 'rotation'
% 'translation'
% 'stretch'
% 'shear'
% 'rotation_translation'
% 'anisotropic_stretch'
% 'vortex'

x=x.'; % transform from N*2 to 2*N


%%
switch mode

    % 1 ROTATION
    case 'rotation'

        omega = pi;
        phi = omega*t;

        A = [cos(phi) -sin(phi);
            sin(phi)  cos(phi)];

        dA = omega*[-sin(phi) -cos(phi);
            cos(phi) -sin(phi)];

        b = [0;0];
        db = [0;0];

        % 2 TRANSLATION
    case 'translation'

        u = [1;0.6];

        A = eye(2);
        dA = zeros(2);

        b = u*t;
        db = u;

        % 3 ISOTROPIC STRETCH
    case 'stretch'

        alpha = 0.5;

        A  = (1+alpha*t)*eye(2);
        dA = alpha*eye(2);

        b  = [0;0];
        db = [0;0];

        % 4 SHEAR
    case 'shear'

        k = 0.8;

        A = [1 k*t;
            0 1];

        dA = [0 k;
            0 0];

        b  = [0;0];
        db = [0;0];

        % 5 ROTATION + TRANSLATION
    case 'rotation_translation'

        omega = pi;
        phi   = omega*t;

        A  = [cos(phi) -sin(phi);
            sin(phi)  cos(phi)];

        dA = omega*[-sin(phi) -cos(phi);
            cos(phi) -sin(phi)];

        u  = [0.8;0.3];

        b  = u*t;
        db = u;

        % 6 ANISOTROPIC STRETCH
    case 'anisotropic_stretch'

        alpha = 0.7;

        S  = [1+alpha*t 0;
            0 1-alpha*t];

        dS = [alpha 0;
            0 -alpha];

        phi = pi/6;

        R = [cos(phi) -sin(phi);
            sin(phi)  cos(phi)];

        A  = R*S*R';
        dA = R*dS*R';

        b  = [0;0];
        db = [0;0];

        % 7 VORTEX (NONLINEAR)
    case 'vortex'

        r     = sqrt(x(1,:).^2 + x(2,:).^2);
        theta = atan2(x(2,:),x(1,:));

        beta  = 1.0;

        theta_new = theta + beta*r*t;

        x_new = [r .* cos(theta_new); r.*sin(theta_new)];

        % velocity
        dtheta_dt = beta*r;

        dx_dt = [-r .* sin(theta_new) .* dtheta_dt; r .* cos(theta_new) .* dtheta_dt];

    case 'inflate'
        % center of inflation
        c = mean(x,2);     % or specify manually

        % inflation rate
        beta = 0.5;

        alpha = 1 + beta*t; % could be something more complex
        dalpha = beta;

        A = alpha*eye(2);
        dA = dalpha*eye(2);
        
        b = c - A*c;
        db = -dA*c;
        % transformed curve
        % x_new = c + alpha*(x - c);
        

        % velocity
        % dx_dt = dalpha*(x - c);


end

%% affine transformations
if ~strcmp(mode,'vortex')

    x_new = A*x + b;

    dx_dt = dA*x + db;


end

x_new = x_new.';
dx_dt = dx_dt.';
