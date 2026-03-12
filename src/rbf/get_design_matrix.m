function A = get_design_matrix(s_samples,epsilon)
Nsamples = length(s_samples);

A = zeros(Nsamples,Nsamples);
if (length(epsilon)==1)
    epsilon = epsilon*ones(Nsamples,1);
end

for i = 1:Nsamples
    for j = 1:Nsamples
        % distance between s_i and s_j along the curve
        gamma_ij = abs(s_samples(i) - s_samples(j)); 
        % For closed curve, use min distance (wrap around 2*pi)
        gamma_ij = min(gamma_ij, 2*pi - gamma_ij);
        % Fill A
        A(i,j) = rbf(gamma_ij,epsilon(j));
    end
end
