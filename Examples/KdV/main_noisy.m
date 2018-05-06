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
%%%    S^3d is applied to analyze the KdV equation with two solitons:   %%%
%%%                    u{t} = -0.000484uu{x} + u{xxx}                   %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
warning off
addpath ../../Functions
%% Generate data

% load kdv.mat
% load Time & Space.mat
data = load('kdv.mat');
u = real(data.u);
dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t;
x = data.x;

% add noise
rng('default');
rng(0);
u = u + 0.01*std(reshape(u,1,size(u,1)*size(u,2)))*randn(size(u,1),size(u,2));

%% Estimate derivatives

% Method: polynomial interpolation
t_th = 1; % order of time derivative
x_th = 3; % maximum order of spatial derivative

parameter.deg = 6;           % degree of polynomial to use 
parameter.x_num_to_fit = 10; % number of points to use in polynomial interpolation for spatial derivatives
parameter.t_num_to_fit = 25; % number of points to use in polynomial interpolation for time derivative

[derivative]= make_input_poly( u,x,x_th,parameter);        
y0 =  make_y_poly( u,t,t_th,parameter );
      
input = derivative.derivative;
input_name = derivative.name;

%% Builds a dictionary matrix

polyorder = 2; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(input(:,1),input(:,2:end),input_name(1),input_name(2:end),polyorder);
        
%% Number of snapshots: 18750

index_all = zeros(size(u,1)-2*parameter.x_num_to_fit,size(u,2)-2*parameter.t_num_to_fit);
for i = 1:size(index_all,1)
    for j = 1:size(index_all,2)
        index_all(i,j) = i+(j-1)*(size(u,1)-2*parameter.x_num_to_fit);
    end
end

index = [];
for i = 1001:1250
    index = [index,index_all(26+floor((i-1)*0.042):100+floor((i-1)*0.042),i)];
end

index = reshape(index,1,size(index,1)*size(index,2));

theta = theta0(index,:);
y = y0(index);

%% Compute Sparse regression: Sparse Bayesian Approach

% normalization
Mreg = zeros(size(theta,2),1);
T = zeros(size(theta,1),size(theta,2));
for i = 1:size(theta,2)
    Mreg(i,1) = 1/max(abs(theta(:,i)));
    T(:,i) = Mreg(i)*theta(:,i);   
end

% dimensionality reduction
[U,S,V] = svd(T);
T = U(:,1:size(theta,2))'*T;
y = U(:,1:size(theta,2))'*y;
clear U;

lambda = 70.8;  % the regularization parameter 
MAXITER = 5;     % maximum number of inner iterations
w = tac_reconstruction(y, T, lambda,MAXITER);
w = w(:,end).*Mreg; % identified parameters corresponding to the basis function in Theta

%% print result

threshold = 1e-6; % filter those cofficients in w being less than the threshold 
y_name = 'u{t}';
fprintf('\n%s = ', y_name);
for i = 1:size(w,1)
    if abs(w(i))<threshold
        w(i) = 0;
    else
        if w(i)<0
            fprintf('%.6f%s', w(i),theta_name{i});
        else
            fprintf('+');
            fprintf('%.6f%s', w(i),theta_name{i});
        end
    end    
end

%% print error

% RMS Error
err = abs([(-0.000484 -w(6))/0.000484,  (-1 - w(7))]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err)*100,std(err)*100);
