function [x_rbf,y_rbf,dxds_rbf,dyds_rbf] = eval_rbf(cx,cy,s,s_samples,epsilon,jac)

Ns    = length(s); % number of evaluation points
x_rbf = zeros(Ns,1);
y_rbf = zeros(Ns,1);


if (nargin>5)
    if (jac==true)
        dxds_rbf = zeros(Ns,1);
        dyds_rbf = zeros(Ns,1);
    end
else
    jac=false;
end


for i = 1:Ns
    gamma_vec = abs(s(i) - s_samples);
    gamma_vec = min(gamma_vec, 2*pi - gamma_vec);  % wrap-around
    if (jac==true)
        [phi,dphi] = rbf(gamma_vec,epsilon,jac);
    else
        [phi,~] = rbf(gamma_vec,epsilon,jac);
    end
    x_rbf(i) = sum(cx .* phi); % sum over all samples
    y_rbf(i) = sum(cy .* phi);
    if (jac==true)
        dxds_rbf(i) = sum(cx .* dphi);
        dyds_rbf(i) = sum(cy .* dphi);
    end
end