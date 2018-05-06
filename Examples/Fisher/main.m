% This file is the main part of S^3d method
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
%%%           S^3d is applied to analyze the Fisher equation:           %%%
%%%                     u{t} = u - u^2 + 0.1u{xx}                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
addpath ../../Functions
%% Generate data

% load fisher.mat
% load Time & Space.mat
data = load('fisher.mat');
u = real(data.u);
dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t;
x = data.x;

%% Estimate derivatives

% Method: finite difference
t_th = 1; % order of time derivative
x_th = 3; % maximum order of spatial derivative
   
[derivative]= make_input_fd( u,dx,x_th);       
y0  = make_y_fd( u,dt,t_th );
        
input = derivative.derivative;
input_name = derivative.name;

%% Builds a dictionary matrix

polyorder = 3; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(input(:,1),input(:,2:end),input_name(1),input_name(2:end),polyorder);
        
%% Number of snapshots: 10000

rng('default');
rng(0);
index = randperm(19000,10000)+11000;
index = reshape(index,1,size(index,1)*size(index,2));

theta = theta0(index,:);
y = y0(index,:);

%% Compute Sparse regression: Sparse Bayesian Approach

% normalization
T = zeros(size(theta,1),size(theta,2));
Mreg = zeros(size(theta,2),1);
for i = 1:size(theta,2)
    Mreg(i,1) = 1/max(abs(theta(:,i)));
    T(:,i) = Mreg(i)*theta(:,i);    
end

lambda = 0.04; % the regularization parameter 
MAXITER = 5;   % maximum number of inner iterations
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
            fprintf('%.4f%s', w(i),theta_name{i});
        else
            fprintf('+');
            fprintf('%.4f%s', w(i),theta_name{i});
        end
    end    
end

%% print error

% RMS Error
err = abs([(1 -w(2)), (-1 -w(3)), (0.1 - w(6))/0.1]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err)*100,std(err)*100);
