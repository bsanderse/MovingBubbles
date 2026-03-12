%% Rotating curve in Cartesian coordinates
% parameterise a curve with x(s,t) and y(s,t), where s is the parameter
% along the curve and t is time. 
% in this code, x and y are treated as separate variables
% in main.m this has been changed to a single variable, which allows
% easier extension to 3D
% this code is also restricted to rotations, whereas main.m allows general
% transforms

clear; clc; close all;
addpath('libs');
addpath('rbf');
addpath('curves');

%% Simulation parameters
s_min     = 0;
s_max     = 2*pi; % parameter along curve
Ns        = 100;   % number of evaluation points
t_min     = 0;
t_max     = 1;
Nt        = 101; % number of time steps
Nsamples  = 20; % number of sample points used to construct RBF
Nt_plot   = 5; % number curves to plot

%% rotation
phi_dot = pi;       % constant rotation speed (rad/s)
phi0 = 0;

%% set up s and t vectors
s = linspace(s_min, s_max, Ns).';   % angles at which to evaluate r
t = linspace(t_min, t_max, Nt).';          % time vector
t_plot = linspace(t_min, t_max, Nt_plot).'; 

% choose sampling points at t=0
% s_samples = s_max*rand(Nsamples,1);         % random angles
s_samples = linspace(0,s_max,Nsamples).';
% s_samples = 0.5*s_max*(legpts(Nsamples)+1); % legendre points on [-1,1]
% s_samples = 0.5*s_max*(chebpts(Nsamples,1)+1);



%% Preallocate location matrices
x = zeros(Ns, Nt);
y = zeros(Ns, Nt);
dx_dt = zeros(Ns, Nt);
dy_dt = zeros(Ns, Nt);

x_samples = zeros(Nsamples, Nt);
y_samples = zeros(Nsamples, Nt);

dx_dt_samples = zeros(Nsamples, Nt);
dy_dt_samples = zeros(Nsamples, Nt);


%% Compute x,y and derivatives 
% at the evaluation points - for plotting purposes
[x(:,1),y(:,1)] = bubble(s);
% at the sampling points - to later update the coefficients when marching
% in time
[x_samples(:,1),y_samples(:,1)] = bubble(s_samples);

IC = [x(:,1) y(:,1)].';
IC_samples = [x_samples(:,1) y_samples(:,1)].';

for k = 1:Nt
    phi = phi0 + phi_dot * t(k); % current rotation angle

    % rotation matrix
    dphi = phi-phi0;
    R    = [cos(dphi) -sin(dphi); sin(dphi) cos(dphi)];

    XYZ    = R*IC;
    x(:,k) = XYZ(1,:);
    y(:,k) = XYZ(2,:);

    Rdot = phi_dot*[-sin(dphi) -cos(dphi); cos(dphi) -sin(dphi)];
    XYZ_dot = Rdot*IC;
    
    dx_dt(:,k) = XYZ_dot(1,:);
    dy_dt(:,k) = XYZ_dot(2,:);

    XYZ_samples    = R*IC_samples;
    x_samples(:,k) = XYZ_samples(1,:);
    y_samples(:,k) = XYZ_samples(2,:);

    XYZ_dot_samples = Rdot*IC_samples;
    
    dx_dt_samples(:,k) = XYZ_dot_samples(1,:);
    dy_dt_samples(:,k) = XYZ_dot_samples(2,:);


    % check time derivative with finite difference
    % if (k>1)
    %     dx_dt_fd(:,k) = (x(:,k) - x(:,k-1)) /(t(k) - t(k-1));
    % end
end

%% Plot x over time for selected s values
figure;
plot(t, x(1,:), 'LineWidth', 1.5); hold on;
plot(t, x(floor(end/2),:), 'LineWidth', 1.5);
plot(t, x(end,:), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('x(s,t)');
legend(sprintf('s = %.2f rad', s(1)), ...
       sprintf('s = %.2f rad', s(floor(end/2))), ...
       sprintf('s = %.2f rad', s(end)));
title('x location of rotating curve at different s');
grid on;

%% Plot dx/dt over time for selected s values
figure;
plot(t, dx_dt(1,:), 'LineWidth', 1.5); hold on;
plot(t, dx_dt(floor(end/2),:), 'LineWidth', 1.5);
plot(t, dx_dt(end,:), 'LineWidth', 1.5);
% plot(t,dr_dt_fd(1,:))
xlabel('Time (s)');
ylabel('dx/dt (s,t)');
legend(sprintf('s = %.2f rad', s(1)), ...
       sprintf('s = %.2f rad', s(floor(end/2))), ...
       sprintf('s = %.2f rad', s(end)));
title('dx/dt at of rotating curve at different s');
grid on;

%% Plot the true curve at different times
figure; hold on; axis equal;
colors = jet(Nt_plot);  % color map for different time instances

for k = 1:Nt_plot   
    [val,loc] = min(abs(t-t_plot(k)));

    plot(x(:,loc), y(:,loc), 'Color', colors(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('t = %.2f s', t(loc)));
end

xlabel('x'); ylabel('y'); axis equal; grid on;
title('Rotating bubble at different time instances');
legend('show');
xlim([-1.5 1.5])
ylim([-1.5 1.5])


%% Define RBF kernel
epsilon = 1;                         % RBF width parameter
% epsilon = mean(abs(diff(s_samples)));

%% Construct design matrix A
A = get_design_matrix(s_samples,epsilon);


%% Solve for coefficients at t=0
cx = A \ x_samples(:,1);  % solve A*c = x
cy = A \ y_samples(:,1);  % solve A*c = y

%% Evaluate RBF reconstruction at t=0
% theta_eval = linspace(0, theta_max, 200);
x_rbf = zeros(Ns,Nt);
y_rbf = zeros(Ns,Nt);

[x_rbf(:,1), y_rbf(:,1)] = eval_rbf(cx,cy,s,s_samples,epsilon);


%% Plot RBF reconstruction and original shape and sampling points

figure; hold on; axis equal; grid on;

plot(x(:,1),y(:,1),'LineWidth', 2,'DisplayName', 'Exact');
plot(x_samples(:,1), y_samples(:,1), 'o', 'MarkerSize', 8, 'DisplayName', 'Sampling points');
plot(x_rbf(:,1), y_rbf(:,1), '-', 'LineWidth', 2, 'DisplayName', 'RBF reconstruction');

xlabel('x'); ylabel('y'); 
title('RBF fit to curve');
legend show;
xlim([-1.5 1.5])
ylim([-1.5 1.5])
box on

%% time dependent coefficients
% we assume we have a fit at t=0, which gives the coefficients cx, cy
% then, for t>0 we assume we know dx/dt at the sampling points, and we use
% this to update c via
% dc/dt = A \ dx/dt
% and subsequently, we find dx/dt at the evaluation points, and then update
% the position with Forward Euler

dx_dt_rbf = zeros(Ns,Nt);
dy_dt_rbf = zeros(Ns,Nt);

% x_rbf     = zeros(Ns,Nt);
% x_samples  = zeros(Nsamples,Nt);
% dx_dt_samples  = zeros(Nsamples,Nt);

% x_rbf(:,1) = bubble(s);
figure; hold on; axis equal; grid on;

for k = 1:Nt

    % phi = phi0 + phi_dot * t(k); % current rotation angle

    % assume we have the time derivatives at the sampling points
    dcx_dt = A\dx_dt_samples(:,k);
    dcy_dt = A\dy_dt_samples(:,k);

    % compute derivatives at evaluation points
    [dx_dt_rbf(:,k), dy_dt_rbf(:,k)] = eval_rbf(dcx_dt,dcy_dt,s,s_samples,epsilon);


    % find new positions of the evaluation points
    if (k>1)
        % forward Euler
        x_rbf(:,k) = x_rbf(:,k-1) + (t(k)-t(k-1))*dx_dt_rbf(:,k-1);
        y_rbf(:,k) = y_rbf(:,k-1) + (t(k)-t(k-1))*dy_dt_rbf(:,k-1);
        
        % plot(x_rbf(:,k),y_rbf(:,k),'-');
    end

    % alternatively, we could compute the new positions of the sample points:
    % x_samples(:,k) = x_samples(:,k-1) + (t(k)-t(k-1))*dx_dt_samples(:,k-1);
    % then compute the coefficients
    % cx = A\x_samples(:,k);
    % and then evaluate the rbf at the evaluation points
    % x_rbf = eval_rbf(cx ,... )

end

%% 

% imid = floor(Nt/2);
plot(x_rbf(:,1),y_rbf(:,1),'-', 'LineWidth', 2, 'DisplayName', 'RBF reconstruction at t=0');
plot(x_rbf(:,end),y_rbf(:,end),'-', 'LineWidth', 2, 'DisplayName', sprintf('RBF reconstruction at t = %.2f', t(end)));
legend
box on
xlim([-1.5 1.5])
ylim([-1.5 1.5])