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
%%%              S^3d is applied to analyze the NLS equation:           %%%
%%%                u{t} = 10/3iu - 10/3iu|u|^2 - 0.3iu{xx}              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
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

%% Estimate derivatives

% Method: Pade
t_th = 1; % order of time derivative
x_th = 3; % maximum order of spatial derivative
       
[derivative_r] = make_input_pade(real(u),dx,x_th);  
[derivative_i] = make_input_pade(imag(u),dx,x_th);  

y0_r = make_y_pade(real(u),dt,t_th);
y0_i = make_y_pade(imag(u),dt,t_th);
y0 = y0_r + y0_i*sqrt(-1);

input = derivative_r.derivative+derivative_i.derivative*sqrt(-1);
input_name = derivative_r.name; 

%% Builds a dictionary matrix

derdata = input(:,2:end);
dername = input_name(2:end);

data = [input(:,1) reshape(abs(u),size(u,1)*size(u,2),1)];
dataname = {'u','|u|'};

polyorder = 3; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(data,derdata,dataname,dername,polyorder);

%% Number of snapshots: 10000

index_all = zeros(size(u,1),size(u,2));

for i = 1:size(index_all,1)
    for j = 1:size(index_all,2)
        index_all(i,j) = i+(j-1)*size(u,1);
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

lambda = 0.001; % the regularization parameter 
MAXITER = 5;  % maximum number of inner iterations
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
