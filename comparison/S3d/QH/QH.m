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
%%%              S^3d is applied to analyze the QH equation:            %%%
%%%                      u{t} = 0.5iu{xx} - i(x^2/2)u                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear 
close all
addpath ../../../Functions

load('data1.mat');
load('data2.mat');
load('data3.mat');

theta1 = [real(theta1) -imag(theta1);imag(theta1) real(theta1) ];
y1 = double([real(y1);imag(y1)]);

theta2 = [real(theta2) -imag(theta2);imag(theta2) real(theta2) ];
y2 = [real(y2);imag(y2)];

theta3 = [real(theta3) -imag(theta3);imag(theta3) real(theta3) ];
y3 = [real(y3);imag(y3)];

%normalization
T1 = zeros(size(theta1,1),size(theta1,2));
Mreg1 = zeros(size(theta1,2),1);
for i = 1:size(theta1,2)
    Mreg1(i,1) = 1/max(abs(theta1(:,i)));
    T1(:,i) = Mreg1(i)*theta1(:,i);    
end
% 
T2 = zeros(size(theta2,1),size(theta2,2));
Mreg2 = zeros(size(theta2,2),1);
for i =  1:size(theta2,2)
    Mreg2(i,1) = 1/max(abs(theta2(:,i)));
    T2(:,i) = Mreg2(i)*theta2(:,i);
end

T3 = zeros(size(theta3,1),size(theta3,2));
Mreg3 = zeros(size(theta3,2),1);
for i =  1:size(theta3,2)
    Mreg3(i,1) = 1/max(abs(theta3(:,i)));
    T3(:,i) = Mreg3(i)*theta3(:,i);
end

MAXITER = 5;    %maximum number of inner iterations
lambda1 = 0.15;  %the regularization parameter 
w1 = tac_reconstruction(y1, T1, lambda1,MAXITER);
w1 = w1(:,end).*Mreg1; %identified parameters corresponding to the basis function in Theta
% 
lambda2 = 0.15;  %the regularization parameter 
w2 = tac_reconstruction(y2, T2, lambda2,MAXITER);
w2 = w2(:,end).*Mreg2; %identified parameters corresponding to the basis function in Theta

lambda3 = 0.16;  %the regularization parameter 
w3 = tac_reconstruction(y3, T3, lambda3,MAXITER);
w3 = w3(:,end).*Mreg3; %identified parameters corresponding to the basis function in Theta

%print result and error

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
err1 = abs([(0.5 -imag(w1(12)))/0.5,(-1 -imag(w1(7)))]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err1)*100,std(err1)*100);
% 
for i = 1:size(w2,1)
    if abs(w2(i))<threshold
        w2(i) = 0;
    end
end

w2 = w2(1:40) + sqrt(-1)*w2(41:80);
fprintf('\ndata2:\n%s = ', y_name);
for i = 1:size(w2,1)
    if abs(w2(i))~=0
        str = strcat('+(',num2str(w2(i),'%.4f'),')');
        fprintf(str);
        fprintf('%s',theta_name{i});
    end
end

err2 = abs([(0.5 -imag(w2(12)))/0.5,(-1 -imag(w2(7)))]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err2)*100,std(err2)*100);

for i = 1:size(w3,1)
    if abs(w3(i))<threshold
        w3(i) = 0;
    end
end

w3 = w3(1:40) + sqrt(-1)*w3(41:80);
fprintf('\ndata3:\n%s = ', y_name);
for i = 1:size(w3,1)
    if abs(w3(i))~=0
        str = strcat('+(',num2str(w3(i),'%.4f'),')');
        fprintf(str);
        fprintf('%s',theta_name{i});
    end
end

err3 = abs([(0.5 -imag(w3(12)))/0.5,(-1 -imag(w3(7)))]);
fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err3)*100,std(err3)*100);

% err = (err1 + err2 + err3)/3;
% fprintf('\nerror: mea = %.4f%%, stad = %.4f%%\n',mean(err)*100,std(err)*100);