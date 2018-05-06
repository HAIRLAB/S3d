function [u_diff] = pade_c4(u,h)

% @author: Ye Yuan
%This function computes the derivative of given data by using the Pade 
%method.
%
%Inputs:
%      u      : one-dimensional data
%      h      : step of sampling
% 
%Output:
%      u_diff : first-order derivative of u
%

n = length(u);

i = 1;
b(i) = 1.0;
c(i) = 2.0;
r(i) = (-5.0*u(i) + 4.0*u(i+1) + u(i+2))/(2.0*h);

for i = 2:n-1
    a(i) = 1.0;
    b(i) = 4.0;
    c(i) = 1.0;
    r(i) = (3.0/h)*(u(i+1)-u(i-1));
end 

i = n;
a(i) = 2.0;
b(i) = 1.0;
r(i) = (-5.0*u(i) + 4.0*u(i-1) + u(i-2))/(-2.0*h);

u_diff = tdma(a,b,c,r,1,n);
u_diff = u_diff.';

end

function [x] = tdma(a,b,c,r,s,e)
%------------------------------------------------------------------!
 
% forward elimination phase
for i=s+1:e
    b(i) = b(i) - a(i)/b(i-1)*c(i-1);
    r(i) = r(i) - a(i)/b(i-1)*r(i-1);
end 
% backward substitution phase 
x(e) = r(e)/b(e);
for i=e-1:-1:s
    x(i) = (r(i)-c(i)*x(i+1))/b(i);
end 

end