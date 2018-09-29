% 
% @author: Ye Yuan
% 
% utilizing Douglas-Rachford algorithm to analyze the NLS equation:            
%                   u{t} = 0.5iu{xx} + iu|u|^2
%

clear
addpath ../../Functions

load('data1.mat');

theta1 = theta01;
y1 = y01;


T1 = zeros(size(theta1,1),size(theta1,2));
Mreg1 = zeros(size(theta1,2),1);
for i = 1:size(theta1,2)
    Mreg1(i,1) = 1/max(abs(theta1(:,i)));
    T1(:,i) = Mreg1(i)*theta1(:,i);
end
%



mu =0.5;
ga = 0.5;

MaxIt = 5000;
lambda = 500;
wall =[];
par =[];
err = [];
num = 1;
start_o = 1e-4;
tol = 1e-7;
for mu = 0.1:0.1:0.9
    for ga = 0.1:0.1:0.9
        for lam_i = 1:12
            for lam_j = 1:9
                lambda = start_o*10^lam_i + start_o*10^lam_i*(lam_j-1);
                
                
                par(num,:) =[ga,lambda,mu];
                [ga,lambda,mu]
                num=num + 1;
                w1w=DouglasRachfordxt(theta1,y1,ga,lambda,mu,MaxIt,tol);
                w1 = w1w(:,end).*Mreg1; %identified parameters corresponding to the basis function in Theta
                for ddd = 1:size(w1,1)
                    if abs(w1(ddd))<1e-6
                        w1(ddd) = 0;
                    end
                end
                
                %w1 = w1(1:40) + sqrt(-1)*w1(41:80);
                wall = [wall w1];
                err1 = mean(abs([(0.5 -imag(w1(12)))/0.5,(1 -imag(w1(9)))]));
                err = [err err1];
            end
        end
    end
end

