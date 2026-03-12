function [x_rbf,dxds_rbf] = eval_rbf_x(c,s,s_samples,epsilon,jac)

Ns    = length(s); % number of evaluation points
x_rbf = zeros(Ns,2,1);


if (nargin>5)
    if (jac==true)
        dxds_rbf = zeros(Ns,1);
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
    x_rbf(i,:) = sum(phi.*c); % sum over all samples
    
    if (jac==true)
        dxds_rbf(i,:) = sum(dphi.*c);
        
    end
end