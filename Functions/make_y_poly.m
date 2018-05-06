function [ut] = make_y_poly(u,t,order,parameter)

% @author: Ye Yuan
%
%This function computes the time derivative on LHS of PDE by using the 
%polynomial interpolation method.
%
%Inputs:
%      u    : two-dimensional data 
%      t    : t-coordinates
%      order: order of time derivative
%      parameter.x_num_to_fit: number of points to use in polynomial 
%                              interpolation for spatial derivatives
%      parameter.t_num_to_fit: number of points to use in polynomial 
%                              interpolation for time derivatives
%      parameter.deg         : degree of polynomial to use 
%Output:
%      ut   : order-th time derivative of u
%

x_num = parameter.x_num_to_fit ;
t_num = parameter.t_num_to_fit ;
deg   = parameter.deg ;

[r,v] = size(u);
ut = zeros(r-2*x_num,v-2*t_num);
num = (r-2*x_num)*(v-2*t_num);

for i = 1: size(ut,1)
    for k = 1:size(ut,2)
        tmp = polyfit(t(k:k+2*t_num)',u(i+x_num,k:k+2*t_num),deg);
        for m = 1:order
            tmp = polyder(tmp );
        end
        ut(i,k) = polyval(tmp,t(k+t_num));
    end
end

ut = reshape(ut,num,1);


