function [derivative] = make_input_fd(varargin)

% @author: Ye Yuan
%This function computes the spatial derivatives in the dictionary functions
%by using the finite difference method.
% 
%Usage 1: [derivative] = make_input_fd(u,dx,order)
%Inputs: 
%      u    : 2-dimensional data
%      dx   : the space step of sampling
%      order: maximum order of spatial derivative 
%Output:
%      derivative.derivative: spatial derivative matrix, including 
%           derivative terms of order-th order, (order-1)-th order... 1st 
%           order, and u
%      derivative.name      : names of the derivative terms
%%%%%
%
%Usage 2: [derivative] = make_input_fd(u,dx,order,dir)
%Inputs:
%      u    : 3-dimensional data
%      h    : the space step of sampling
%      order: order of spatial derivative 
%      dir  : name of the spatial variable, should be 'x' or 'y'
%Output:
%      derivative: order-th spatial derivative of u with respect to dir
%

if nargin == 3
    
    u     = varargin{1};
    dx    = varargin{2};
    order = varargin{3};
    
    ux = zeros(size(u));
    num = size(ux,1)*size(ux,2);

    derivative.derivative = reshape(u,num,1);
    name_num = 1;
    name = {'u'};
    namex = {'x','xx','xxx'};

    for j = 1:order
        for i = 1: size(u,2)
            ux(:,i) = finitediff(u(:,i),dx,j);
        end
        derivative.derivative = [derivative.derivative  reshape(ux,num,1)];
        name_tmp = ['u{' namex{j} '}'];
        name_num = name_num+1;
        name{name_num} = name_tmp;    
    end
    derivative.name = name;

elseif nargin == 4
    
    u     = varargin{1};
    h     = varargin{2};
    order = varargin{3};
    dir   = varargin{4};
    
    derivative = zeros(size(u));    
    if strcmp(dir,'x')
        for j = 1:size(u,2)
            for k = 1:size(u,3)
                derivative(:,j,k) = finitediff(u(:,j,k),h,order);
            end
        end
    
    elseif strcmp(dir,'y')
        for i = 1:size(u,1)
            for k = 1:size(u,3)
                derivative(i,:,k) = finitediff(u(i,:,k),h,order);
            end
        end
    end
    
end
