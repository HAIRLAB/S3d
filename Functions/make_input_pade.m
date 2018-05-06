function [derivative] = make_input_pade(varargin)

% @author: Ye Yuan
%This function computes the spatial derivatives in the dictionary functions
%by using Pade method.
% 
%Usage 1: [derivative] = make_input_pade(u,h,order)
%Inputs: 
%      u    : 2-dimensional data
%      h    : the space step of sampling
%      order: maximum order of spatial derivative 
%Output:
%      derivative.derivative: spatial derivative matrix, including 
%           derivative terms of order-th order, (order-1)-th order... 1st 
%           order, and u
%      derivative.name      : names of the derivative terms
%%%%%
%
%Usage 2: [derivative] = make_input_pade(u,dx,order,dir)
%Inputs:
%      u    : 3-dimensional data
%      h   : the space step of sampling
%      order: order of spatial derivative 
%      dir  : name of the spatial variable, should be 'x' or 'y'
%Output:
%      derivative: order-th spatial derivative of u with respect to dir
%

if nargin == 3
    
    u     = varargin{1};
    h     = varargin{2};
    order = varargin{3};
    
    name = {'u'};
    namex = {'x','xx','xxx','xxxx','xxxxx'};

    for i = 1:order
        name_tmp = ['u{' namex{i} '}']; 
        name{i+1} = name_tmp;
    end
    derivative.name = name;
    tmp = u;
    derivative.derivative = zeros(size(u,1)*size(u,2),order+1);
    derivative.derivative(:,1) = reshape(u,size(u,1)*size(u,2),1); 
    for i = 1:order
        for j = 1:size(u,2)
            tmp(:,j) = pade_c4(tmp(:,j),h);
        end
        derivative.derivative(:,i+1) = reshape(tmp,size(tmp,1)*size(tmp,2),1);
    end    

elseif nargin ==4
    
    u     = varargin{1};
    h     = varargin{2};
    order = varargin{3};
    dir   = varargin{4};
    
    derivative = zeros(size(u));
    if strcmp(dir,'x')
        derivative = u;
        for m = 1:order
            for j = 1:size(u,2)
                for k = 1:size(u,3)
                    derivative(:,j,k) = pade_c4(derivative(:,j,k),h);
                end
            end
        end

    elseif strcmp(dir,'y')
        derivative = u;
        for m = 1:order
            for i = 1:size(u,1)
                for k = 1:size(u,3)
                    derivative(i,:,k) = pade_c4(derivative(i,:,k),h);
                end
            end
        end
    end   
    
end
        