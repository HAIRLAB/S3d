function [up]=pade_c4(u,h)

up = zeros(size(u));
n = length(u);

i=1;
b(i) = 1.0;
c(i) = 2.0;
r(i) = (-5.0*u(i) + 4.0*u(i+1) + u(i+2))/(2.0*h);

for i=2:n-1
a(i) = 1.0;
b(i) = 4.0;
c(i) = 1.0;
r(i) = (3.0/h)*(u(i+1)-u(i-1));
end 

i=n;
a(i) = 2.0;
b(i) = 1.0;
r(i) = (-5.0*u(i) + 4.0*u(i-1) + u(i-2))/(-2.0*h);

up=tdma(a,b,c,r,1,n);
up=up';
end