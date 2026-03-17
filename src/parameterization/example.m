% Example: noisy sphere
N = 500;
theta = rand(N,1)*pi;
phi = rand(N,1)*2*pi;

X = [sin(theta).*cos(phi), ...
     sin(theta).*sin(phi), ...
     cos(theta)];

% Run Isomap
[Y, D] = isomap(X, 10);

% coordinates s, q are Y(:,1) and Y(:,2)
% in other words, for a given point with index i and coordinate X(i,),
% the corresponding parameters are Y(i,1) and Y(i,2)

% Plot
figure;
scatter(Y(:,1), Y(:,2), 20, theta, 'filled');
title('Isomap embedding (s,q)');
axis equal;

