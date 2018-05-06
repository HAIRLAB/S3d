function [derivative] = make_input_poly(u,x,order,parameter)

% @author: Ye Yuan
%
%This function computes the spatial derivatives in the dictionary functions 
%by using the polynomial interpolation method.
%
%Inputs:
%      u    : two-dimensional data 
%      x    : x-coordinates
%      order: maximum order of spatial derivative
%      parameter.x_num_to_fit: number of points to use in polynomial 
%                              interpolation for spatial derivatives
%      parameter.t_num_to_fit: number of points to use in polynomial 
%                              interpolation for time derivatives
%      parameter.deg         : degree of polynomial to use 
%
%Output:
%      derivative.derivaive  : spatial derivative matrix, including 
%           derivative terms of order-th order, (order-1)-th order... 1st 
%           order, and u
%      derivative.name       : names of the derivative terms

deg = parameter.deg ;
x_num = parameter.x_num_to_fit ;
t_num = parameter.t_num_to_fit ;

[r,v] = size(u);
ux = zeros(r-2*x_num,v-2*t_num);
num = (r-2*x_num)*(v-2*t_num);
derivative.derivative = reshape(u(1+x_num:end-x_num,1+t_num:end-t_num),num,1);

name_num = 1;
name = {'u'};
namex = {'x','xx','xxx','xxxx','xxxxx'};

Fit = zeros(size(ux,1)*size(ux,2),deg+1);
for i = 1: size(ux,2)
    for k = 1:size(ux,1)
        Fit((i-1)*size(ux,1)+k,:) = polyfit(x(k:k+2*x_num)',u(k:k+2*x_num,i+t_num),deg);
    end
end

for j = 1:order
    for i = 1: size(ux,2)
        for k = 1:size(ux,1)            
            tmp = Fit((i-1)*size(ux,1)+k,:);
            for m = 1:j
                tmp = polyder(tmp);
            end
            ux(k,i) = polyval(tmp,x(k+x_num));
        end
    end
    derivative.derivative = [derivative.derivative  reshape(ux,num,1)];
    name_tmp1 = ['u{' namex{j} '}' ];
    name_num = name_num+1;
    name{name_num} = name_tmp1;
end

derivative.name = name;
