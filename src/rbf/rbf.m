function [phi,dphidx] = rbf(x,epsilon,jac)
% epsilon can be scalar or vector, in the latter case it should have the
% dimension of x
phi =  exp(-(x./epsilon).^2);  % Gaussian RBF
if (nargin>=3)
    if (jac==true)
        dphidx = -2*x./(epsilon.^2).*phi;
    else
        dphidx = 0;
    end
end
