%% Transforming a curve in Cartesian coordinates
% parameterise a curve with x(s,t) and y(s,t), where s is the parameter
% along the curve and t is time.
% apply a linear or nonlinear rotation to it
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
Nt        = 41; % number of time steps
Nsamples  = 50; % number of sample points used to construct RBF
Nt_plot   = 5; % number curves to plot

curve_type = 'bubble'; % 'bubble', 'ellipse', ..

%% transformation
% 'rotation'
% 'translation'
% 'stretch'
% 'shear'
% 'rotation_translation'
% 'anisotropic_stretch'
% 'vortex' 
% 'inflate'

transform_type = 'anisotropic_stretch';

%% set up s and t vectors
s = linspace(s_min, s_max, Ns).';   % parameter locations at which to evaluate x
t = linspace(t_min, t_max, Nt).';          % time vector
t_plot = linspace(t_min, t_max, Nt_plot).';

% choose sampling points at t=0
% s_samples = s_max*rand(Nsamples,1);         % random angles
s_samples = linspace(0,s_max,Nsamples).';     % uniform samples
% s_samples = 0.5*s_max*(legpts(Nsamples)+1); % legendre points on [-1,1]
% s_samples = 0.5*s_max*(chebpts(Nsamples,1)+1);



%% Preallocate location matrices
x = zeros(Ns, 2, Nt); % contains both x and y locations
dx_dt = zeros(Ns, 2, Nt);

x_samples = zeros(Nsamples, 2, Nt);

dx_dt_samples = zeros(Nsamples, 2, Nt);


%% Compute x,y and derivatives
% at the evaluation points - for plotting purposes:
x(:,:,1) = parametric_curve(s,curve_type);

% at the sampling points - to later update the coefficients when marching
% in time:
x_samples(:,:,1) = parametric_curve(s_samples,curve_type);

for k = 1:Nt

    % apply transform to initial condition
    [x(:,:,k), dx_dt(:,:,k)] = transform(x(:,:,1),t(k),transform_type);

    [x_samples(:,:,k), dx_dt_samples(:,:,k)] = transform(x_samples(:,:,1),t(k),transform_type);

    % check time derivative with finite difference
    % if (k>1)
    %     dx_dt_fd(:,k) = (x(:,k) - x(:,k-1)) /(t(k) - t(k-1));
    % end
end


%% Plot x over time for selected s values
figure;
plot(t, squeeze(x(1,1,:)), 'LineWidth', 1.5); hold on;
plot(t, squeeze(x(floor(end/2),1,:)), 'LineWidth', 1.5);
plot(t, squeeze(x(end,1,:)), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('x(s,t)');
legend(sprintf('s = %.2f rad', s(1)), ...
    sprintf('s = %.2f rad', s(floor(end/2))), ...
    sprintf('s = %.2f rad', s(end)));
title('x location of transformed curve at different s');
grid on;

%% Plot dx/dt over time for selected s values
figure;
plot(t, squeeze(dx_dt(1,1,:)), 'LineWidth', 1.5); hold on;
plot(t, squeeze(dx_dt(floor(end/2),1,:)), 'LineWidth', 1.5);
plot(t, squeeze(dx_dt(end,1,:)), 'LineWidth', 1.5);
% plot(t,dr_dt_fd(1,:))
xlabel('Time (s)');
ylabel('dx/dt (s,t)');
legend(sprintf('s = %.2f rad', s(1)), ...
    sprintf('s = %.2f rad', s(floor(end/2))), ...
    sprintf('s = %.2f rad', s(end)));
title('dx/dt at of transformed curve at different s');
grid on;

%% Plot the true curve at different times
figure; hold on; axis equal;
colors = jet(Nt_plot);  % color map for different time instances

for k = 1:Nt_plot
    [val,loc] = min(abs(t-t_plot(k)));

    plot(x(:,1,loc), x(:,2,loc), 'Color', colors(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('t = %.2f s', t(loc)));
end

xlabel('x'); ylabel('y'); axis equal; grid on;
title('Transformed parametric curve at different time instances');
legend('show');
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])


%% Define RBF kernel
% epsilon = 1;                         % RBF width parameter
% arc length: 
% epsilon = mean(abs(diff(s_samples)));  % heuristic
% epsilon = mean(norm(diff(x_samples),2));  % heuristic

% heuristic based on distance between subsequent points along curve
d = sqrt(sum(diff(x_samples(:,:,1)).^2,2));
h = mean(d);
epsilon = 1.5*h;

% local h and corresponding epsilon
% D = pdist2(x_samples(:,:,1),x_samples(:,:,1));           % pairwise distance
% D(1:Nsamples+1:end) = inf;        % ignore self-distance
% h = min(D,[],2);            % nearest-neighbor spacing
% epsilon = h;

%% Construct design matrix A
A = get_design_matrix(s_samples,epsilon);


%% Solve for coefficients at t=0
c = A \ x_samples(:,:,1);  % solve A*c = x


%% Evaluate RBF reconstruction at t=0
% theta_eval = linspace(0, theta_max, 200);
x_rbf = zeros(Ns,2,Nt);

x_rbf(:,:,1) = eval_rbf_x(c,s,s_samples,epsilon);


%% Plot RBF reconstruction and original shape and sampling points

figure; hold on; axis equal; grid on;

plot(x(:,1,1),x(:,2,1),'LineWidth', 2,'DisplayName', 'Exact');
plot(x_samples(:,1,1), x_samples(:,2,1), 'o', 'MarkerSize', 8, 'DisplayName', 'Sampling points');
plot(x_rbf(:,1,1), x_rbf(:,2,1), '-', 'LineWidth', 2, 'DisplayName', 'RBF reconstruction');

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

dx_dt_rbf = zeros(Ns,2,Nt);
x_samples_num = zeros(Nsamples,2,Nt);
x_samples_num(:,:,1) = x_samples(:,:,1);


for k = 1:Nt

    % assume we have the time derivatives at the sampling points
    dc_dt = A\dx_dt_samples(:,:,k);

    % compute derivatives at evaluation points
    dx_dt_rbf(:,:,k) = eval_rbf_x(dc_dt,s,s_samples,epsilon);

    % find new positions of the evaluation points
    if (k>1)
        % forward Euler for now - can be easily be improved
        x_rbf(:,:,k) = x_rbf(:,:,k-1) + (t(k)-t(k-1))*dx_dt_rbf(:,:,k-1);

        % note we can also compute the new positions of the sample points:
        x_samples_num(:,:,k) = x_samples_num(:,:,k-1) + (t(k)-t(k-1))*dx_dt_samples(:,:,k-1);
        % and then we can compute the coefficients
        % cx = A\x_samples(:,k);
        % and then evaluate the rbf at the evaluation points
        % x_rbf = eval_rbf(cx ,... )
    end

end

%%
figure; hold on; axis equal; grid on;

tmid = floor(Nt/2)+1;
plot(x_samples(:,1,1), x_samples(:,2,1), 'o', 'LineWidth', 1,'MarkerSize', 8, 'DisplayName', 'exact sample points at t=0');
plot(x_rbf(:,1,1),x_rbf(:,2,1),'-', 'LineWidth', 2, 'DisplayName', 'RBF reconstruction at t=0');

plot(x_samples(:,1,tmid),x_samples(:,2,tmid),'o', 'LineWidth', 1, 'DisplayName', sprintf('exact sample points at t = %.2f', t(tmid)));
plot(x_samples_num(:,1,tmid),x_samples_num(:,2,tmid),'x', 'LineWidth', 1, 'DisplayName', sprintf('approximated sample points at t = %.2f', t(tmid)));
plot(x_rbf(:,1,tmid),x_rbf(:,2,tmid),'-', 'LineWidth', 2, 'DisplayName', sprintf('RBF reconstruction at t = %.2f', t(tmid)));

plot(x_samples(:,1,end),x_samples(:,2,end),'o', 'LineWidth', 1, 'DisplayName', sprintf('exact sample points at t = %.2f', t(end)));
plot(x_samples_num(:,1,end),x_samples_num(:,2,end),'x', 'LineWidth', 1, 'DisplayName', sprintf('approximated sample points at t = %.2f', t(end)));
plot(x_rbf(:,1,end),x_rbf(:,2,end),'-', 'LineWidth', 2, 'DisplayName', sprintf('RBF reconstruction at t = %.2f', t(end)));


legend('Location','northwest');
box on
xlim([-2.5 2])
ylim([-1.75 2.5])
xlabel('x')
ylabel('y')
title('Bubble deforming under anisotropic stretch')