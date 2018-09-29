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
%%%     Ye Yuan, "Discovery of Partial Differential Equations from      %%%
%%%                        Spatiotemporal Data"                         %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              S^3d is applied to analyze the NLS equation:           %%%
%%%                      u{t} = 0.5iu{xx} + iu|u|^2                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
warning off
addpath ../../../Functions
%% Generate data

% load nls.mat
% load Time & Space.mat
data = load('nls.mat');
u = data.usol;
dt = data.t(2) - data.t(1);
dx = data.x(2) - data.x(1);
t = data.t;
x = data.x;

%add noise
rng('default');
rng(0);
u = u + 0.01/sqrt(2)*std(reshape(real(u),1,size(u,1)*size(u,2)))*randn(size(u))+0.01/sqrt(2)*std(reshape(imag(u),1,size(u,1)*size(u,2)))*randn(size(u));

%% Estimate derivatives

%Method: polynomial interpolation
t_th = 1; %order of time derivative
x_th = 3; %maximum order of spatial derivative

parameter.deg = 6;           %degree of polynomial to use 
parameter.x_num_to_fit = 20; %number of points to use in polynomial interpolation for spatial derivatives
parameter.t_num_to_fit = 10; %number of points to use in polynomial interpolation for time derivative

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

%% Number of snapshots: 

index_all = zeros(size(u,1)-2*x_num,size(u,2)-2*t_num);
for i = 1:size(index_all,1)
    for j = 1:size(index_all,2)
        index_all(i,j) = i+(j-1)*(size(u,1)-2*x_num);
    end
end

% index = index_all(101:400,111:150);
% index = reshape(index,1,size(index,1)*size(index,2));
%%
index = index_all(131:250,121:150);
index = index_all(191:250,371:400);
index = reshape(index,1,size(index,1)*size(index,2));

theta01 = theta0(index,:);
y01 = y0(index,:);
% save('data1','theta01','y01','theta_name');

%%
% save('data1','theta1','y1','theta_name');
% theta00 = [real(theta0) -imag(theta0);imag(theta0) real(theta0) ];
% y00 = [real(y0);imag(y0)];
theta1 = [real(theta01) -imag(theta01);imag(theta01) real(theta01) ];
y1 = [real(y01);imag(y01)];
%normalization
T1 = zeros(size(theta1,1),size(theta1,2));
Mreg1 = zeros(size(theta1,2),1);
for i = 1:size(theta1,2)
    Mreg1(i,1) = 1/max(abs(theta1(:,i)));
    T1(:,i) = Mreg1(i)*theta1(:,i);    
end

%dimensionality reduction
[U,~,~] = svd(T1);
T1 = U(:,1:80)'*T1;
y1 = U(:,1:80)'*y1;

clear U;
%
%%
MAXITER =10;    %maximum number of inner iterations
lambda1 =999.9;  %the regularization parameter 
w1w = tac_reconstruction(y1, T1, lambda1,MAXITER);
w1 = w1w(:,end).*Mreg1; %identified parameters corresponding to the basis function in Theta
threshold = 1e-6; %filter those cofficients in w being less than the threshold
y_name = 'u{t}';

for i = 1:size(w1,1)
    if abs(w1(i))<threshold
        w1(i) = 0;
    end
end

w1 = w1(1:40) + sqrt(-1)*w1(41:80);
fprintf('\ndata1:\n%s = ', y_name);
for i = 1:size(w1,1)
    if abs(w1(i))~=0
        str = strcat('+(',num2str(w1(i),'%.4f'),')');
        fprintf(str);
        fprintf('%s',theta_name{i});
    end
end

%RMS Error
err1 = abs([(0.5 -imag(w1(12)))/0.5,(1 -imag(w1(9)))]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err1)*100,std(err1)*100);

