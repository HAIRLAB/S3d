% 
% @author: Ye Yuan
% 
% utilizing Douglas-Rachford algorithm to analyze the QH equation:            
%                   u{t} = 0.5iu{xx} - i(x^2/2)u
%

clear
clc 
close all

load('data3.mat');

theta1 = [real(theta3) -imag(theta3);imag(theta3) real(theta3) ];
y1 = [real(y3);imag(y3)];

%normalization
T1 = zeros(size(theta1,1),size(theta1,2));
Mreg1 = zeros(size(theta1,2),1);
for i = 1:size(theta1,2)
    Mreg1(i,1) = 1/max(abs(theta1(:,i)));
    T1(:,i) = Mreg1(i)*theta1(:,i);
end

tau =0;
mu = 0.5;
MaxIt = 1e5;
tol = 0;
sigma = 10;
wall =[];
par =[];
err = [];
num = 1;
start_o = 1e-4;
for mu_k = 1:8
    mu = 0.1+mu_k*0.1
    for sigma_i = 1:6
        for sigma_j = 1:9
            sigma = start_o*10^sigma_i + start_o*10^sigma_i*(sigma_j-1);
            for tau_i = 1:6
                for tau_j = 1:9
                    tau = start_o*10^tau_i + start_o*10^tau_i*(tau_j-1);

                    par(num,:) = [mu,sigma,tau];
                    [mu,sigma,tau]
                    num=num + 1;
                    w1w=DouglasRachford(T1,y1,sigma,tau,mu,MaxIt,tol);
                    w1 = w1w(:,end).*Mreg1; %identified parameters corresponding to the basis function in Theta
                    for ddd = 1:size(w1,1)
                        if abs(w1(ddd))<1e-6
                            w1(ddd) = 0;
                        end
                    end
                    
                    w1 = w1(1:40) + sqrt(-1)*w1(41:80);
                    wall = [wall w1];
                    err1 = mean(abs([(0.5 -imag(w1(12)))/0.5,(-1 -imag(w1(7)))]));
                    err = [err err1];
                end
            end
        end
    end
end

%% plot qh3.fig

% origin_data = load('data3');
% theta1 = origin_data.theta3;
% y1 = origin_data.y3;
% 
% theta1 = [real(theta1) -imag(theta1);imag(theta1) real(theta1) ];
% y1 = double([real(y1);imag(y1)]);
% for i = 1:size(wall,2)
%     tmp_w = abs(wall(:,i));
%     nonzeros(i) = length(find(tmp_w~=0));
%     w1 = [real(wall(:,i));imag(wall(:,i))];
%     y_err(i) = norm(y1-theta1*w1,2)^2/norm(y1,2)^2;
%     
% end
% x_index = 1:size(wall,2);
% [AX,H1,H2]=plotyy(x_index,nonzeros,x_index,y_err);
% set(H1,'Marker','.','LineStyle','none',...
% 'Color',[0.0784313753247261 0.168627455830574 0.549019634723663]);
% 
% set(H2,'LineWidth',1.5,...
% 'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
% set(AX(1),'FontSize',20,'LineWidth',1.5,'YColor',[0 0.447 0.741],'YTick',...
% [0 20 40]);
% set(AX(2),'FontSize',20,'Color','none','HitTest','off','YAxisLocation','right','YColor',...
% [0.85 0.325 0.098],'YTick',[0 0.5 1]);
% 
% set(get(AX(1),'Ylabel'),'String','Number of terms') 
% set(get(AX(2),'Ylabel'),'String','Fitting error') 
% xlabel({'Parameters'});