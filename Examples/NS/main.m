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
%%%         S^3d is applied to analyze the Navier Stokes equation:      %%%
%%%              w{t} = 0.01w{xx} + 0.01w{yy} -uw{x} - vw{y}            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
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

num = size(u,1)*size(u,2)*size(u,3);

%% Estimate derivatives

% Method: finite difference
y0 = make_y_fd(w,dt,1);

wx = make_input_fd(w,dx,1,'x');
wy = make_input_fd(w,dy,1,'y');
wxx = make_input_fd(w,dx,2,'x');
wxy = make_input_fd(wx,dy,1,'y');
wyy = make_input_fd(w,dy,2,'y');

%% Builds a dictionary matrix

data = [reshape(u,num,1) reshape(v,num,1) reshape(w,num,1)];
dataname = {'u','v','w'};

derdata = [reshape(wx,num,1) reshape(wy,num,1) reshape(wxx,num,1) reshape(wxy,num,1) reshape(wyy,num,1)];
dername = {'w{x}','w{y}','w{xx}','w{xy}','w{yy}'};

polyorder = 2; % maximum power of polynomial function to be included in Theta
[theta0,theta_name] = build_Theta(data,derdata,dataname,dername,polyorder);

%% Number of snapshots: 10000

xindex = 41:60;
yindex = 21:40;
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

lambda = 0.00008; % the regularization parameter 
MAXITER = 5;      % maximum number of inner iterations
ww = tac_reconstruction(y,T,lambda,MAXITER);
w = ww(:,end).*Mreg; % identified parameters corresponding to the basis function in Theta

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
