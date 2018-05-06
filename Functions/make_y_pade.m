function [ut] = make_y_pade(u,h,order)

% @author: Ye Yuan
%
%This function computes the time derivative on LHS of PDE by using Pade 
%method.
%
%Inputs:
%      u    : data(2 or 3 dimensional)
%      h    : the time step of sampling
%      order: order of time derivative
%
%Output:
%      ut   : order-th time derivative of u
%

if length(size(u)) == 2
    
    ut = u;
    for j = 1:order
        for i = 1:size(u,1)
            ut(i,:) = pade_c4(ut(i,:),h);
        end
    end

    ut = reshape(ut,size(ut,1)*size(ut,2),1);
    
elseif length(size(u)) ==3
    
    ut = u;
    for m = 1:order
        for i = 1:size(u,1)
            for j = 1:size(u,2)
                ut(i,j,:) = pade_c4(ut(i,j,:),h);
            end
        end
    end
    
    ut = reshape(ut,size(ut,1)*size(ut,2)*size(ut,3),1);
    
end

