function [up]=c4(u,h)
%implicit none
%integer:n,i
%step: h
%dimension (1:n): u
up = zeros(size(u));
n = length(u);
for i=3:n-2
up(i) = (-u(i+2)+8.0*u(i+1)-8.0*u(i-1)+u(i-2))/(12.0*h);
end 

% sided difference for i=0 (2nd-order)
i=1;
up(i) = (-3.0*u(i) + 4.0d0*u(i+1) - u(i+2))/(2.0*h);

% sided difference for i=1 (3rd-order)
i=2;
up(i) = (-2.0*u(i-1) - 3.0*u(i) + 6.0*u(i+1) - u(i+2))/(6.0*h);

% sided difference for i=n (2nd-order)
i=n;
up(i) = (-3.0*u(i) + 4.0*u(i-1) - u(i-2))/(-2.0*h);

% sided difference for i=n-1 (3rd-order)
i=n-1;
up(i) = (-2.0*u(i+1) - 3.0*u(i) + 6.0*u(i-1) - u(i-2))/(-6.0*h);

end
