function [x]=tdma(a,b,c,r,s,e)
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