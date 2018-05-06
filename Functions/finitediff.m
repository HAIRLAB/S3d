function [u_diff] = finitediff(u,d,order)

% @author: Ye Yuan
%This function computes the derivative of given data by using the finite 
%difference method.
%
%Inputs:
%      u      : one-dimensional data
%      d      : step of sampling
%      order  : order of derivative
% 
%Output:
%      u_diff : order-th derivative of u
%

u_diff = zeros(size(u));

n = length(u);

switch order
    case 1
        for i = 2 : n-1
            u_diff(i) = (u(i+1) - u(i-1))/(2*d);
        end
        u_diff(1) = (-3.0/2*u(1) + 2*u(2) - u(3)/2) / d;
        u_diff(n) = (3.0/2*u(n) - 2*u(n-1) + u(n-2)/2) / d;
        
    case 2
        for i = 2 : n-1
            u_diff(i) = (u(i+1)-2*u(i)+u(i-1)) / d^2;
        end
        u_diff(1) = (2*u(1) - 5*u(2) + 4*u(3) - u(4)) / d^2;
        u_diff(n) = (2*u(n) - 5*u(n-1) + 4*u(n-2) - u(n-3)) / d^2;
        
    case 3
        for i = 3 : n-2
            u_diff(i) = (u(i+2)/2-u(i+1)+u(i-1)-u(i-2)/2) / d^3;
        end
        u_diff(1) = (-2.5*u(1)+9*u(2)-12*u(3)+7*u(4)-1.5*u(5)) / d^3;
        u_diff(2) = (-2.5*u(2)+9*u(3)-12*u(4)+7*u(5)-1.5*u(6)) / d^3;
        u_diff(n-1) = (2.5*u(n)-9*u(n-1)+12*u(n-2)-7*u(n-3)+1.5*u(n-4)) / d^3;
        u_diff(n) = (2.5*u(n-1)-9*u(n-2)+12*u(n-3)-7*u(n-4)+1.5*u(n-5)) / d^3;
        
end


