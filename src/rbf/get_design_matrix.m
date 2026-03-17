function A = get_design_matrix(s_samples,epsilon,type)
Nsamples = length(s_samples);

% A = zeros(Nsamples,Nsamples);
if (isscalar(epsilon)) % make a vector of epsilon
    epsilon = epsilon*ones(Nsamples,1);
end

%non-vectorized
% for i = 1:Nsamples
%     for j = 1:Nsamples
%         % distance between s_i and s_j along the curve
%         gamma_ij = abs(s_samples(i) - s_samples(j)); 
%         % For closed curve, use min distance (wrap around 2*pi)
%         gamma_ij = min(gamma_ij, 2*pi - gamma_ij);
%         % Fill A
%         A(i,j) = rbf(gamma_ij,epsilon(j),false,type);
%     end
% end

% vectorized
eps_row = epsilon(:).';             % 1 x Nsamples
s = s_samples(:);                   % ensure column
G = abs(s - s.');                   % [Nsamples x Nsamples] pairwise abs diff
G = min(G, 2*pi - G);               % wrap-around distance
A = rbf(G,eps_row,false,type);