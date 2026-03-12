function [x,y]=bubble(s,R,alpha,beta)

if (nargin==1)
    R=1; % radius
end
if (nargin<=2)
    alpha=0.75;
end
if (nargin<=3)
    beta=-0.5;
end


x = R*sin(s);
y = R*(cos(s) + alpha*cos(2*s) + beta*cos(3*s));

