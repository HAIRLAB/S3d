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
%%%   S^3d is applied to analyze the complex Ginzburg-Landau equation   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
addpath ../../Functions
%% Generate data

% load ca06212cb.mat
dataR = load('car06212cb');
dataL = load('cal06212cb');

AR = dataR.car06212cb(:,1) - 1i*dataR.car06212cb(:,2);
AL = dataL.cal06212cb(:,1) + 1i*dataL.cal06212cb(:,2);

dt = 1.5625;   %time step of sampling
dx = 0.458167; %space step of sampling

x_len = 180;
t_len = size(AR,1)/x_len;

AR = reshape(AR,x_len,t_len);
AL = reshape(AL,x_len,t_len);

%% Estimate derivatives

AL_dx = zeros(size(AL));
AL_dxx = zeros(size(AL));

for i = 1:t_len
    AL_dx(:,i) = pade_c4(AL(:,i),dx);
    AL_dxx(:,i) = pade_c4(AL_dx(:,i),dx);    
end

AL_dt = zeros(size(AL));

for i = 1:x_len
    AL_dt(i,:) = pade_c4(AL(i,:),dt);
end

%% Builds a dictionary matrix

AL_dx = reshape(AL_dx,x_len*t_len,1);
AL_dxx = reshape(AL_dxx,x_len*t_len,1);
AL_dt = reshape(AL_dt,x_len*t_len,1);
AR = reshape(AR,x_len*t_len,1);
AL = reshape(AL,x_len*t_len,1);

theta0 = [AL_dx,AL,AL_dxx,abs(AL).^2.*AL,abs(AR).^2.*AL];
y0 = AL_dt;

%% Number of snapshots: 136800

index_all = zeros(x_len,t_len);
for i = 1:size(index_all,1)
    for j = 1:size(index_all,2)
        index_all(i,j) = i+(j-1)*x_len;
    end
end

index = index_all(:,21:780);
index = reshape(index,size(index,1)*size(index,2),1);

theta = theta0(index,:);
y = y0(index,:);

rng('default');
rng(0);
index = randperm(136800,136800);
theta = theta(index,:);
y = y(index);

%4/5 for training
theta_train = theta(1:109440,:);
y_train = y(1:109440);
theta_train = [real(theta_train) -imag(theta_train(:,2:5));imag(theta_train) real(theta_train(:,2:5))];
y_train = [real(y_train);imag(y_train)];

%1/5 for validation
theta_val = theta(109441:end,:);
y_val = y(109441:end);
theta_val = [real(theta_val) -imag(theta_val(:,2:5));imag(theta_val) real(theta_val(:,2:5))];
y_val = [real(y_val);imag(y_val)];

%% Compute Sparse regression: Sparse Bayesian Approach

%normalization
Mreg = zeros(size(theta_train,2),1);
T = zeros(size(theta_train,1),size(theta_train,2));
for i = 1:size(theta_train,2)
    Mreg(i,1) = 1/max(abs(theta_train(:,i)));
    T(:,i) = Mreg(i)*theta_train(:,i);   
end

%dimensionality reduction
stepsize = length(y_train)/8;
y2 = [];
T2 = [];
for i = 1:8
    [U,S,V] = svd(T((i-1)*stepsize+1:i*stepsize,:));
    T2 = [T2;U(:,1:9)'*T((i-1)*stepsize+1:i*stepsize,:)];
    y2 = [y2;U(:,1:9)'*y_train((i-1)*stepsize+1:i*stepsize)];
end
clear U;

T = T2;
y_train = y2;

lambda = zeros(126,1);

for i = 1:7
    for j = 1:18
        lambda((i-1)*18+j) = 0.5*(20-j)*10^(-5-i); %the regularization parameter 
    end
end

threshold = 5e-4; %filter those cofficients in w being less than the threshold
MAXITER = 5;      %maximum number of inner iterations
error = zeros(126,1);
w_end = zeros(9,126);

for i = 1:length(lambda)
    
    w= tac_reconstruction(y_train, T, lambda(i),MAXITER);
    w = w(:,end).*Mreg; %identified parameters corresponding to the basis function in Theta
    
    for j = 1:size(w,1)
        if abs(w(j))<threshold
            w(j) = 0;
        end
    end
    w_end(:,i) = w; 
    error(i) = norm(y_val-theta_val*w,2)^2/norm(y_val,2)^2; %compute fitting error on the validation set
    
end

[minerror,loc] = min(error);

%% obtain w corresponding to the selected lambda

lambda_final = lambda(loc); 

theta = [real(theta) -imag(theta(:,2:5));imag(theta) real(theta(:,2:5))];
y = [real(y);imag(y)];

%normalization
Mreg = zeros(size(theta,2),1);
T = zeros(size(theta,1),size(theta,2));
for i = 1:size(theta,2)
    Mreg(i,1) = 1/max(abs(theta(:,i)));
    T(:,i) = Mreg(i)*theta(:,i);   
end

%dimensionality reduction
stepsize = length(y)/10;
y2 = [];
T2 = [];
for i = 1:10
    [U,S,V] = svd(T((i-1)*stepsize+1:i*stepsize,:));
    T2 = [T2;U(:,1:9)'*T((i-1)*stepsize+1:i*stepsize,:)];
    y2 = [y2;U(:,1:9)'*y((i-1)*stepsize+1:i*stepsize)];
end
clear U;

T = T2;
y = y2;

MAXITER = 5;
w = tac_reconstruction(y, T, lambda_final,MAXITER);
w = w(:,end).*Mreg;

%filter those cofficients in w being less than the threshold
threshold = 5e-4;
for k = 1:size(w,1)
    if abs(w(k))<threshold
        w(k) = 0;
    end
end

%print result
theta_name = {'s','epsilon/tau0','xi0*xi0/tau0','g/tau0','h/tau0','epsilon*c0/tau0','xi0*xi0*c1/tau0','g*c2/tau0','h*c3/tau0'};
for i = 1:9
    fprintf(theta_name{i});
    fprintf(' = %.4f\n',w(i));
end