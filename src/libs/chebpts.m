function x = chebpts(n,kind)

k = (1:n).';

if (kind==1)

    % First kind (T_n) zeros
    x = cos((2*k-1)/(2*n)*pi);
elseif (kind==2)
    % Second kind (U_n) zeros
    x = cos(k/(n+1)*pi);

end
x = sort(x);


end