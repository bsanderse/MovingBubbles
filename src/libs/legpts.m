function x = legpts(n)
% stable numerical computation using eigenvalues of Jacobi matrix

beta = (1:n-1) ./ sqrt((2*(1:n-1)+1).*(2*(1:n-1)-1)); % recurrence betas
J = diag(beta,1) + diag(beta,-1);    % symmetric tridiagonal Jacobi matrix
x = sort(eig(J));   

end