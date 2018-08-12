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
%%%         S^3d is applied to analyze the Navier Stokes equation:      %%%
%%%              w{t} = 0.01w{xx} + 0.01w{yy} -uw{x} - vw{y}            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
warning off
addpath ../../Functions
%% Generate data

% load ns.mat
% load Time & Space.mat
data = load('ns.mat');
dt = data.dt;
dx = data.dx;
dy = data.dy;
u = data.uuu;
v = data.vvv;
w = data.www;

% add noise
rng('default');
rng(0);

num = size(u,1)*size(u,2)*size(u,3);

u = u + 0.01*std(reshape(u,1,num))*randn(size(u));
v = v + 0.01*std(reshape(v,1,num))*randn(size(v));
w = w + 0.01*std(reshape(w,1,num))*randn(size(w));

% use POD to denoise
u = reshape(u,100*100,1001);
[phi_u,~,nbasis_u] = POD(u,99.99);
u = phi_u(:,1:nbasis_u)*phi_u(:,1:nbasis_u)'*u;
clear phi_u;
u = reshape(u,100,100,1001);

v = reshape(v,100*100,1001);
[phi_v,~,nbasis_v] = POD(v,99.99);
v = phi_v(:,1:nbasis_v)*phi_v(:,1:nbasis_v)'*v;
clear phi_v;
v = reshape(v,100,100,1001);

w = reshape(w,100*100,1001);
[phi_w,~,nbasis_w] = POD(w,99.99);
w = phi_w(:,1:nbasis_w)*phi_w(:,1:nbasis_w)'*w;
clear phi_w;
w = reshape(w,100,100,1001);

%% Estimate derivatives

% Method: Pade
y0 = make_y_pade(w,dt,1);

wx  = make_input_pade(w,dx,1,'x');
wy  = make_input_pade(w,dy,1,'y');
wxx = make_input_pade(w,dx,2,'x');
wxy = make_input_pade(wx,dy,1,'y');
wyy = make_input_pade(w,dy,2,'y');

%% Builds a dictionary matrix

derdata = [reshape(wx,num,1) reshape(wy,num,1) reshape(wxx,num,1) reshape(wxy,num,1) reshape(wyy,num,1)];
dername = {'w{x}','w{y}','w{xx}','w{xy}','w{yy}'};

data = [reshape(u,num,1) reshape(v,num,1) reshape(w,num,1)];
dataname = {'u','v','w'};

polyorder = 2; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(data,derdata,dataname,dername,polyorder);

%% Number of snapshots: 10000

xindex = 41:60;
yindex = 41:60;
tindex = 101:125;

theta = zeros(20*20*25,60);
y = zeros(20*20*25,1);

for i = 1:20
    for j = 1:20
        for k = 1:25
            theta(i+(j-1)*20+(k-1)*400,:) = theta0(xindex(i)+(yindex(j)-1)*100+(tindex(k)-1)*10000,:);
            y(i+(j-1)*20+(k-1)*400) = y0(xindex(i)+(yindex(j)-1)*100+(tindex(k)-1)*10000);
        end
    end
end

%% Compute Sparse regression: Sparse Bayesian Approach

% normalization
T = zeros(size(theta,1),size(theta,2));
Mreg = zeros(size(theta,2),1);
for i = 1:size(theta,2)
    Mreg(i,1) = 1/max(abs(theta(:,i)));
    T(:,i) = Mreg(i)*theta(:,i);    
end

% dimensionality reduction
% [U,S,V] = svd(T);
% T = U(:,1:60)'*T;
% y = U(:,1:60)'*y;
% clear U;

lambda = 0.0005; % the regularization parameter 
MAXITER = 10;    % maximum number of inner iterations
w = tac_reconstruction(y,T,lambda,MAXITER);
w = w(:,end).*Mreg;

%% print result

threshold = 1e-6; % filter those cofficients in w being less than the threshold
y_name = 'w{t}';
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
err = abs([(0.01 -w(13))/0.01,(0.01 -w(15))/0.01,(-1 -w(16))  (-1 -w(22))]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err)*100,std(err)*100);
