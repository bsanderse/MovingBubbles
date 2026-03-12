%% Rotating ellipse in polar coordinates
% r(theta,t) = radius at angle theta as ellipse rotates
clear; clc; close all;
addpath('libs');
addpath('rbf');
addpath('curves');

%% Ellipse parameters
a = 2;              % semi-major axis
b = 1;              % semi-minor axis
phi_dot = pi;     % constant rotation speed (rad/s)
phi0 = 0;

%% Simulation parameters
theta_max = pi;
Ntheta = 50; % number of evaluation points
Nt  = 41; % number of time steps
Nsamples = 40; % number of samples used to construct RBF

theta = linspace(0, theta_max, Ntheta).';   % angles at which to evaluate r
t = linspace(0, 1, Nt).';          % time vector

% at t=0
% theta_samples = theta_max*rand(Nsamples,1);         % random angles
% theta_samples = linspace(0,theta_max,Nsamples).';
% theta_samples = 0.5*theta_max*(legpts(Nsamples)+1); % legendre points on [-1,1]
theta_samples = 0.5*theta_max*(chebpts(Nsamples,1)+1);

%% Preallocate radius matrix
r = zeros(Ntheta, Nt);
dr_dt = zeros(Ntheta, Nt);

%% Compute r(theta,t) and dr/dt
for k = 1:Nt
    phi = phi0 + phi_dot * t(k); % current rotation angle

    [r(:,k),dr_dt(:,k)] = ellipse(theta-phi,phi_dot,a,b);
    
    % check time derivative with finite difference
    % if (k>1)
    %     dr_dt_fd(:,k) = (r(:,k) - r(:,k-1)) /(t(k) - t(k-1));
    % end
end

%% Plot radius over time for selected angles
figure;
plot(t, r(1,:), 'LineWidth', 1.5); hold on;
plot(t, r(floor(end/2),:), 'LineWidth', 1.5);
plot(t, r(end,:), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Radius r(\theta,t)');
legend(sprintf('theta = %.2f rad', theta(1)), ...
       sprintf('theta = %.2f rad', theta(floor(end/2))), ...
       sprintf('theta = %.2f rad', theta(end)));
title('Radius of Rotating Ellipse at Different Angles');
grid on;

%% Optional: Plot dr/dt over time
figure;
plot(t, dr_dt(1,:), 'LineWidth', 1.5); hold on;
plot(t, dr_dt(floor(end/2),:), 'LineWidth', 1.5);
plot(t, dr_dt(end,:), 'LineWidth', 1.5);
% plot(t,dr_dt_fd(1,:))
xlabel('Time (s)');
ylabel('dr/dt (\theta,t)');
legend(sprintf('theta = %.2f rad', theta(1)), ...
       sprintf('theta = %.2f rad', theta(floor(end/2))), ...
       sprintf('theta = %.2f rad', theta(end)));
title('Rate of Change of Radius dr/dt at Different Angles');
grid on;

%% Plot the actual ellipse at different times
figure; hold on; axis equal;
colors = jet(Nt);  % color map for different time instances

for k = 1:Nt
    
    % Polar to Cartesian coordinates
    [x,y] = pol2cart(theta,r(:,k));    
   
    plot(x, y, 'Color', colors(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('t = %.2f s', t(k)));
end

xlabel('x'); ylabel('y'); axis equal; grid on;
title('Rotating Ellipse at Different Time Instances');
legend('show');


%% Generate points on the ellipse

[r_samples, ~] = ellipse(theta_samples-phi0,0,a,b);
% r_rand = 1 ./ sqrt((cos(theta_rand - phi0).^2)/a^2 + (sin(theta_rand - phi0).^2)/b^2);

% Cartesian coordinates (optional, just for plotting)
[x_samples,y_samples] = pol2cart(theta_samples,r_samples);

figure; hold on; axis equal; grid on;
plot(x_samples, y_samples, 'o', 'MarkerSize', 8, 'DisplayName', 'Random points');

%% Define RBF kernel
% epsilon = 0.1;                         % RBF width parameter
epsilon = mean(abs(diff(theta_samples)));

%% Construct design matrix A
A = zeros(Nsamples,Nsamples);
for i = 1:Nsamples
    for j = 1:Nsamples
        % angular distance between theta_i and theta_j on the circle
        gamma_ij = abs(theta_samples(i) - theta_samples(j)); 
        % For circular domain, use min distance (wrap around 2*pi)
        gamma_ij = min(gamma_ij, 2*pi - gamma_ij);
        % Fill A
        A(i,j) = rbf(gamma_ij,epsilon);
    end
end

%% Solve for coefficients
c = A \ r_samples;  % solve A*c = r

%% Evaluate RBF reconstruction
% theta_eval = linspace(0, theta_max, 200);
r_rbf = zeros(size(theta));

for k = 1:Ntheta
    gamma_vec = abs(theta(k) - theta_samples);
    gamma_vec = min(gamma_vec, 2*pi - gamma_vec);  % wrap-around
    r_rbf(k) = sum(c .* rbf(gamma_vec,epsilon));
    
end

%% Plot RBF reconstruction
[x_rbf,y_rbf] = pol2cart(theta,r_rbf);

% x_rbf = r_rbf .* cos(theta);
% y_rbf = r_rbf .* sin(theta);
plot(x_rbf, y_rbf, '-', 'LineWidth', 2, 'DisplayName', 'RBF reconstruction');
% hold on
% r_ex = 1 ./ sqrt((cos(theta - phi0).^2)/a^2 + (sin(theta - phi0).^2)/b^2);
r_ex = ellipse(theta-phi0,0,a,b);
[x_ex,y_ex] = pol2cart(theta,r_ex);

plot(x_ex,y_ex,'LineWidth', 2,'DisplayName', 'Exact')
xlabel('x'); ylabel('y'); axis equal;
title('RBF Fit to an Ellipse');
legend show;


%% now, try RBF fit based on the function u~1/r^2
u_rand = 1 ./ (r_samples.^2);


%% design matrix
A_u = A;
% A_u = zeros(Nsamples,Nsamples);
% 
% for i=1:Nsamples
%     for j=1:Nsamples
% 
%         % angular distance between theta_i and theta_j on the circle
%         gamma_ij = abs(theta_rand(i) - theta_rand(j)); 
%         % For circular domain, use min distance (wrap around 2*pi)
%         gamma_ij = min(gamma_ij, 2*pi - gamma_ij);
%         % Fill A
%         A_u(i,j) = rbf(gamma_ij,epsilon);
% 
%     end
% end

%% solve for coefficients
c = A_u \ u_rand;

%% evaluation grid
% theta = linspace(0,2*pi,400)';

u_rbf = zeros(Ntheta,1);

for i=1:Ntheta
    
    phi = zeros(Nsamples,1);
    
    for j=1:Nsamples
        
        gamma_ij = abs(theta(i) - theta_samples(j));
        % For circular domain, use min distance (wrap around 2*pi)
        gamma_ij = min(gamma_ij, 2*pi - gamma_ij);

        phi(j) = rbf(gamma_ij,epsilon); %exp(-(dtheta/epsilon)^2);
        
    end
    
    u_rbf(i) = phi' * c;
    
end

%% recover radius
r_rbf = 1 ./ sqrt(u_rbf);

%% convert to Cartesian for plotting
[x_rbf,y_rbf] = pol2cart(theta,r_rbf);


%% plot
plot(x_rbf,y_rbf,'r--','LineWidth',2,'DisplayName','RBF on 1/r^2')


%% time dependent coefficients
% we assume we have a fit at t=0;
% then, for t>0 we assume we know dr/dt at the sampling points, and we use
% this to update c via
% dc/dt = A \ dr/dt
% and subsequently, we find dr/dt at the evaluation points, and then update
% the position with Forward Euler

figure; hold on; axis equal; grid on;

dr_dt_rbf = zeros(Ntheta,Nt);
r_rbf = zeros(Ntheta,Nt);
r_samples  = zeros(Nsamples,Nt);
dr_dt_samples  = zeros(Nsamples,Nt);

r_rbf(:,1) = ellipse(theta-phi0,0,a,b);


for k = 1:Nt

    phi = phi0 + phi_dot * t(k); % current rotation angle

    [~,dr_dt_samples(:,k)] = ellipse(theta_samples-phi,phi_dot,a,b);

    dc_dt = A\dr_dt_samples(:,k);

    % evaluate at Ntheta angles
    for i = 1:Ntheta
        gamma_vec = abs(theta(i) - theta_samples);
        gamma_vec = min(gamma_vec, 2*pi - gamma_vec);  % wrap-around
        dr_dt_rbf(i,k) = sum(dc_dt .* rbf(gamma_vec,epsilon));
    end
    % plot(theta,dr_dt_rbf(:,k));

    % find new positions
    if (k>1)
        % forward Euler
        r_rbf(:,k) = r_rbf(:,k-1) + (t(k)-t(k-1))*dr_dt_rbf(:,k-1);
        [x_rbf,y_rbf] = pol2cart(theta,r_rbf);
        plot(x_rbf,y_rbf,'-');
    end

end


%% we could repeat the time-dependent simulation for the fit on u instead of r (to be done)
