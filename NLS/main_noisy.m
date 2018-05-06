% This file is the main part of S^3d method(with 1% noise)
%
% Copyright (c) 2018 Ye Yuan
% All rights reserved.
%
% Code by Ye Yuan
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%         This Matlab File uses the S^3 method described in           %%%
%%%    Ye Yuan, "Machine Discovery of Partial Differential Equations    %%%
%%%                        from Spatiotemporal Data"                    %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              S^3d is applied to analyze the NLS equation:           %%%
%%%                u{t} = 10/3iu - 10/3iu|u|^2 - 0.3iu{xx}              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
warning off
addpath ../../Functions
%% Generate data

% load nls.mat
% load Time & Space.mat
data = load('nls.mat');
u = data.usol;
dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t;
x = data.x;

% add noise
rng('default');
rng(0);
u = u + 0.01/sqrt(2)*std(reshape(real(u),1,size(u,1)*size(u,2)))*randn(size(u))+0.01/sqrt(2)*std(reshape(imag(u),1,size(u,1)*size(u,2)))*randn(size(u));

%% Estimate derivatives

% Method: polynomial interpolation
t_th = 1; % order of time derivative
x_th = 3; % maximum order of spatial derivative

parameter.deg = 5;           % degree of polynomial to use 
parameter.x_num_to_fit = 25; % number of points to use in polynomial interpolation for spatial derivatives
parameter.t_num_to_fit = 10; % number of points to use in polynomial interpolation for time derivative

[derivative_r]= make_input_poly( real(u),x,x_th,parameter);
[derivative_i]= make_input_poly( imag(u),x,x_th,parameter);

y0_r =  make_y_poly( real(u),t,t_th,parameter );     
y0_i =  make_y_poly( imag(u),t,t_th,parameter );  
y0 = y0_r + y0_i*sqrt(-1);    

input = derivative_r.derivative+derivative_i.derivative*sqrt(-1);
input_name = derivative_r.name;

%% Builds a dictionary matrix

x_num = parameter.x_num_to_fit;
t_num = parameter.t_num_to_fit;

derdata = input(:,2:end);
dername = input_name(2:end);

data = [input(:,1),reshape(abs(u(1+x_num:end-x_num,1+t_num:end-t_num)),(size(u,1)-2*x_num)*(size(u,2)-2*t_num),1)];
dataname = {'u','|u|'};

polyorder = 3; %maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(data,derdata,dataname,dername,polyorder);

%% Number of snapshots: 10000

index_all = zeros(size(u,1)-2*x_num,size(u,2)-2*t_num);
for i = 1:size(index_all,1)
    for j = 1:size(index_all,2)
        index_all(i,j) = i+(j-1)*(size(u,1)-2*x_num);
    end
end
index = index_all(151:250,151:250);
index = reshape(index,1,size(index,1)*size(index,2));

theta = theta0(index,:);
y = y0(index,:);
theta = [real(theta) -imag(theta);imag(theta) real(theta) ];
y = [real(y);imag(y)];

%% Compute Sparse regression: Sparse Bayesian Approach

% normalization
T = zeros(size(theta,1),size(theta,2));
Mreg = zeros(size(theta,2),1);
for i = 1:size(theta,2)
    Mreg(i,1) = 1/max(abs(theta(:,i)));
    T(:,i) = Mreg(i)*theta(:,i);    
end

% dimensionality reduction
[U,S,V] = svd(T);
T = U(:,1:size(T,2))'*T;
y = U(:,1:size(T,2))'*y;
clear U;

lambda = 0.5; % the regularization parameter 
MAXITER = 5; % maximum number of inner iterations
w = tac_reconstruction(y, T, lambda,MAXITER);
w = w(:,end).*Mreg; % identified parameters corresponding to the basis function in Theta

%% print result

threshold = 1e-6; % filter those cofficients in w being less than the threshold

for i = 1:size(w,1)
    if abs(w(i))<threshold
        w(i) = 0;
    end
end

w = w(1:40) + sqrt(-1)*w(41:80);
y_name = 'u{t}';
fprintf('\n%s = ', y_name);
for i = 1:size(w,1)
    if abs(w(i))~=0
        str = strcat('+(',num2str(w(i),'%.4f'),')');
        fprintf(str);
        fprintf('%s',theta_name{i});
    end
end

%% print error

% RMS Error
err = abs([(-10/3 -imag(w(2)))/(10/3),(10/3 - imag(w(9)))/(10/3),(0.3 - imag(w(12)))/0.3]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err)*100,std(err)*100);
