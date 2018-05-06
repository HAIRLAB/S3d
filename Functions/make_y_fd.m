function [ut] = make_y_fd(u,dt,order)

% @author: Ye Yuan
%
%This function computes the time derivative on LHS of PDE by using the 
%finite difference method.
%
%Inputs:
%      u    : data(2 or 3 dimensional)
%      dt   : the time step of sampling
%      order: order of time derivative
%
%Output:
%      ut   : order-th time derivative of u
%

if length(size(u)) == 2
    
    ut = zeros(size(u));
    num = size(ut,1)*size(ut,2);
    for i = 1:size(u,1)
        ut(i,:) = finitediff( u(i,:),dt,order);
    end
    ut = reshape(ut,num,1);
    
elseif length(size(u)) == 3
    
    ut = zeros(size(u));
    num = size(ut,1)*size(ut,2)*size(ut,3);  
    for i = 1:size(u,1)
        for j = 1:size(u,2)
            ut(i,j,:) = finitediff(u(i,j,:),dt,order);
        end
    end
    ut = reshape(ut,num,1);
    
end


