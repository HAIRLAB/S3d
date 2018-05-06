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
%%%            S^3d is applied to analyze the KS equation:              %%%
%%%                   u{t} = -uu{x} - u{xx} - u{xxxx}                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
warning off;
addpath ../../Functions
%% Generate data

% load ks.mat
% load Time & Space.mat
data = load('ks.mat');
u = real(data.uu);
dt = data.tt(2) - data.tt(1);
dx = data.x(2) - data.x(1);
t = data.tt;
x = data.x;

% add noise
rng('default');
rng(0);
u = u + 0.01*std(reshape(u,1,size(u,1)*size(u,2)))*randn(size(u,1),size(u,2));

% use POD to denoise
[phi,~,nbasis] = POD(u,99.99);
u = phi(:,1:nbasis)*phi(:,1:nbasis)'*u;

%% Estimate derivatives

% Method: polynomial interpolation
t_th = 1; % order of time derivative
x_th = 5; % maximum order of spatial derivative

parameter.deg = 8;           % degree of polynomial to use 
parameter.x_num_to_fit = 45; % number of points to use in polynomial interpolation for spatial derivatives
parameter.t_num_to_fit = 40; % number of points to use in polynomial interpolation for time derivative

[derivative] = make_input_poly( u,x,x_th,parameter);       
y0 = make_y_poly( u , t, t_th, parameter );
               
input = derivative.derivative;
input_name = derivative.name;

%% Builds a dictionary matrix

polyorder = 5; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(input(:,1),input(:,2:end),input_name(1),input_name(2:end),polyorder);

%% Number of snapshots: 59210

index_all = zeros(size(u,1)-2*parameter.x_num_to_fit,size(u,2)-2*parameter.t_num_to_fit);
for i = 1:size(index_all,1)
    for j = 1:size(index_all,2)
        index_all(i,j) = i+(j-1)*(size(u,1)-2*parameter.x_num_to_fit);
    end
end

index = index_all(1:955,end-61:end);

theta = theta0(index,:);
y = y0(index,:);

%% Compute Sparse regression: Sparse Bayesian Approach

% normalization
T=zeros(size(theta,1),size(theta,2));
Mreg = zeros(size(theta,2),1);
for i = 1:size(theta,2)
    Mreg(i,1) = 1/max(abs(theta(:,i)));
    T(:,i) = Mreg(i)*theta(:,i);    
end

% dimensionality reduction
T1 = T(1:size(T,1)/2,:);
y1 = y(1:size(T,1)/2);
[U,~,~] = svd(T1);
T1 = U(:,1:size(theta,2))'*T1;
y1 = U(:,1:size(theta,2))'*y1;

T2 = T(size(T,1)/2+1:end,:);
y2 = y(size(T,1)/2+1:end);
[U,~,~] = svd(T2);
T2 = U(:,1:size(theta,2))'*T2;
y2 = U(:,1:size(theta,2))'*y2;

clear U;

T = [T1;T2];
y = [y1;y2];

lambda = 185; % the regularization parameter 
MAXITER = 10; % maximum number of inner iterations
w= tac_reconstruction(y, T, lambda,MAXITER);
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
            fprintf('%.4f%s', w(i),theta_name{i});
        else
            fprintf('+');
            fprintf('%.4f%s', w(i),theta_name{i});
        end
    end    
end

%% print error

% RMS Error
err = abs([(-1 -w(8)),  (-1 - w(10)), (-1 - w(12))]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err)*100,std(err)*100);
