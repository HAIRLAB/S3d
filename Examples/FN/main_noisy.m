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
%%%            S^3d is applied to analyze the FN equation:              %%%
%%%            u{t} = -0.2u + 1.2u^2 - u^3 + u{xx} - w  ...(1)          %%%
%%%            w{t} = 0.002u - 0.001w                   ...(2)          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
warning off;
addpath ../../Functions

%% Generate data

% load fn.mat
% load Time & Space.mat
data = load('fn.mat');
u = data.u;
w = data.w;
t = data.t;
x = data.x;
dt = data.t(2)-data.t(1);
dx = data.x(2)-data.x(1);

% add noise
rng('default');
rng(0);
u = u + 0.01*std(reshape(u,1,size(u,1)*size(u,2)))*randn(size(u,1),size(u,2));
w = w + 0.01*std(reshape(w,1,size(w,1)*size(w,2)))*randn(size(w,1),size(w,2));

%% Esitimate derivatives

% Method: polynomial interpolation
t_th = 1;   % order of time derivative of u and w
x_th_u = 3; % maximum order of spatial derivative of u
x_th_w = 1; % maximum order of spatial derivative of w

parameter1.deg = 6; % degree of polynomial to use 
parameter1.x_num_to_fit = 30; % number of points to use in polynomial interpolation for spatial derivatives
parameter1.t_num_to_fit = 15; % number of points to use in polynomial interpolation for time derivative of u

[derivative_u]= make_input_poly( u,x,x_th_u,parameter1);
[derivative_w]= make_input_poly( w,x,x_th_w,parameter1);

y0_u = make_y_poly( u ,t,t_th,parameter1 );

parameter2.deg = 4; % degree of polynomial to use 
parameter2.x_num_to_fit = 30;
parameter2.t_num_to_fit = 30; % number of points to use in polynomial interpolation for time derivative of w

y0_w = make_y_poly( w,t,t_th,parameter2);

input_u = derivative_u.derivative;
input_name_u = derivative_u.name;

input_w = derivative_w.derivative;
input_name_w = {'w','w{x}'};

%% Builds a dictionary matrix

polyorder_u = 4; % maximum power of polynomial function in terms of u to be included in Theta 
[theta0_u,theta_name_u] = build_Theta(input_u(:,1),input_u(:,2:end),input_name_u(1),input_name_u(2:end),polyorder_u);
polyorder_w = 2; % maximum power of polynomial function in terms of w to be included in Theta 
[theta0_w,theta_name_w] = build_Theta(input_w(:,1),input_w(:,2:end),input_name_w(1),input_name_w(2:end),polyorder_w);

theta0 = [theta0_u,theta0_w(:,2:end)];
theta_name = [theta_name_u,theta_name_w(2:end)];

%% Number of snapshots: 18000/27000

index_all = zeros(size(u,1)-2*parameter1.x_num_to_fit,size(u,2)-2*parameter1.t_num_to_fit);
for i = 1:size(index_all,1)
    for j = 1:size(index_all,2)
        index_all(i,j) = i+(j-1)*(size(u,1)-2*parameter1.x_num_to_fit);
    end
end

index_u = index_all(181:280,1:180);
index_u = reshape(index_u,1,size(index_u,1)*size(index_u,2));

index_w = index_all(151:300,1:180);
index_w = reshape(index_w,1,size(index_w,1)*size(index_w,2));

y_u = y0_u(index_u,:);
theta_u = theta0(index_u,:);

y_w = y0_w(index_w,:);
index_w = index_w + (512-30*2)*15;
theta_w = theta0(index_w,:);

%% Compute Sparse regression: Sparse Bayesian Approach

% normalization
T_u = zeros(size(theta_u,1),size(theta_u,2));
Mreg_u = zeros(size(theta_u,2),1);
for i = 1:size(theta_u,2)
    Mreg_u(i,1) = 1/max(abs(theta_u(:,i)));
    T_u(:,i) = Mreg_u(i)*theta_u(:,i);    
end

T_w = zeros(size(theta_w,1),size(theta_w,2));
Mreg_w = zeros(size(theta_w,2),1);
for i = 1:size(theta_w,2)
    Mreg_w(i,1) = 1/max(abs(theta_w(:,i)));
    T_w(:,i) = Mreg_w(i)*theta_w(:,i);    
end

% dimensionality reduction
[U,~,~] = svd(T_u);
T_u = U(:,1:25)'*T_u;
y_u = U(:,1:25)'*y_u;

[U,~,~] = svd(T_w);
T_w = U(:,1:25)'*T_w;
y_w = U(:,1:25)'*y_w;
clear U;

MAXITER = 5;      % maximum number of inner iterations
lambda_u = 0.019; % the regularization parameter
w_u = tac_reconstruction(y_u, T_u, lambda_u,MAXITER);
w_u = w_u(:,end).*Mreg_u; % identified parameters corresponding to the basis function in Theta

lambda_w = 6e-5; 
w_w = tac_reconstruction(y_w, T_w, lambda_w,MAXITER);
w_w = w_w(:,end).*Mreg_w;

%% print result

threshold = 1e-6; % filter those cofficients in w being less than the threshold
y_name_u = 'u{t}';
fprintf('\n%s = ', y_name_u);
for i = 1:size(w_u,1)
    if abs(w_u(i))<threshold
        w_u(i) = 0;
    else
        if w_u(i)<0
            fprintf('%.4f%s', w_u(i),theta_name{i});
        else
            fprintf('+');
            fprintf('%.4f%s', w_u(i),theta_name{i});
        end
    end    
end

y_name_w = 'w{t}';
fprintf('\n%s = ', y_name_w);
for i = 1:size(w_w,1)
    if abs(w_w(i))<threshold
        w_w(i) = 0;
    else
         if w_w(i)<0
            fprintf('%.6f%s', w_w(i),theta_name{i});
        else
            fprintf('+');
            fprintf('%.6f%s', w_w(i),theta_name{i});
        end
    end    
end

%% print error

% RMS Error
err = abs([(-0.2 -w_u(2))/0.2,  (1 - w_u(7)),(1.2 - w_u(3))/1.2,(-1 - w_u(4)),(-1 - w_u(21)),(0.002 -w_w(2))/0.002, (-0.001 - w_w(21))/0.001]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err)*100,std(err)*100);
