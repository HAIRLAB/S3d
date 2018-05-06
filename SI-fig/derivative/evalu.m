%% SI-Derivative: Estimate the derivative of analytical solution
%                            to KdV by different method 


%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file verify the performance of different method used to estimate 
%%%     the derivative presented in
%%% Ye Yuan, "Supplementary Information: Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%%% Method: Based on the analytical function
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used functions: 
%         finitediff: 2-order difference scheme
%         c4: 4-order difference scheme
%         pade_c4: a 4th-order compact Pade scheme
%         singlesd: analytical solution 
%         singlesd1: 1-order analytically derivative 
%         singlesd2: 2-order analytically derivative 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the data
c1 = 5.0;
c2 = 1.0;
n = 50;
m = 50;
x = linspace(-10, 12,n);
dt = 0.025;
dx = x(2)-x(1);
t = linspace(dt,m*dt,m);

U0 = zeros(n,m);
U1 = zeros(n,m);
U2 = zeros(n,m);
U1d1 = zeros(n,m);
U1d2 = zeros(n,m);
U1d3 = zeros(n,m);
for i=1:n
    for j=1:m
        U0(i,j) = singlesoliton(x(i),t(j),c1,-3);
        U1(i,j) = singlesoliton(x(i),t(j),c1,-3);
        U2(i,j) = singlesoliton(x(i),t(j),c2,-1);
        U1d1(i,j)=singlesd(x(i),t(j),c1,-3);
        U1d2(i,j)=singlesd1(x(i),t(j),c1,-3);
        U1d3(i,j)=singlesd2(x(i),t(j),c1,-3);
    end
end

%Test1
U1=U1(:,1);
U1d1=U1d1(:,1);
uf=finitediff(U1,dx,1);
u4=c4(U1,dx);
up=pade_c4(U1,dx);


%Test2
U1d2=U1d2(:,1);
uf2=finitediff(uf,dx,1);
u42=c4(u4,dx);
up2=pade_c4(up,dx);

%Test3
U1d3=U1d3(:,1);
uf3=finitediff(uf2,dx,1);
u43=c4(u42,dx);
up3=pade_c4(up2,dx);

%Show results
%%
subplot(1,3,1)
%clf; hold on;
plot(U1d1,'o','MarkerSize',9,'LineWidth',1.5,'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
hold on;
plot(uf,'square','MarkerSize',8,'LineWidth',3,'Color',[0.388235300779343 0.545098066329956 0.854901969432831]);
hold on;
plot(u4,'.','MarkerSize',25,'LineWidth',3,'Color' ,[0 0 0]);
hold on;
plot(up,'^','MarkerSize',8,'Color' ,[0.0784313753247261 0.168627455830574 0.549019634723663]);
%legend('Analytic','2-order','4-order','pade-4','Location','NorthWest','Location','NorthEast');
%xlabel('x')
%ylabel('1-order derivative')
title([ 'u_x']);
set(gca,'linewidth',3);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20);

%figure(2); clf; hold on;
subplot(1,3,2)
%clf; hold on;
plot(U1d2,'o','MarkerSize',9,'LineWidth',1.5,'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
hold on;
plot(uf2,'square','MarkerSize',8,'LineWidth',3,'Color',[0.388235300779343 0.545098066329956 0.854901969432831]);
hold on;
plot(u42,'.','MarkerSize',25,'LineWidth',3,'Color' ,[0 0 0]);
hold on;
plot(up2,'^','MarkerSize',8,'Color' ,[0.0784313753247261 0.168627455830574 0.549019634723663]);
%legend('Analytic','2-order','4-order','pade-4','Location','NorthWest','Location','NorthEast');
%xlabel('x')
%ylabel('2-order derivative')
title([ 'u_{xx}']);
set(gca,'linewidth',3);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20);



%figure(3); clf; hold on;
subplot(1,3,3)
plot(U1d3,'o','MarkerSize',9,'LineWidth',1.5,'Color',[0.635294139385223 0.0784313753247261 0.184313729405403]);
hold on;
plot(uf3,'square','MarkerSize',8,'LineWidth',3,'Color',[0.388235300779343 0.545098066329956 0.854901969432831]);
hold on;
plot(u43,'.','MarkerSize',25,'LineWidth',3,'Color' ,[0 0 0]);
hold on;
plot(up3,'^','MarkerSize',8,'Color' ,[0.0784313753247261 0.168627455830574 0.549019634723663]);
legend('Analytic','2-order','4-order','pade-4','Location','NorthWest','Location','NorthEast');
%xlabel('x')
%ylabel('3-order derivative')
title([ 'u_{xxx}']);
set(gca,'linewidth',3);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20);


%%
figure
subplot(1,3,1)
plot(uf-U1d1,'LineWidth',3,'Color',[0.388235300779343 0.545098066329956 0.854901969432831]);
hold on;
plot(u4-U1d1,'LineWidth',3,'Color' ,[0 0 0]);
hold on;
plot(up-U1d1,'LineWidth',3,'Color' ,[0.0784313753247261 0.168627455830574 0.549019634723663]);
set(gca,'linewidth',3);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20);



subplot(1,3,2)
plot(uf2-U1d2,'LineWidth',3,'Color',[0.388235300779343 0.545098066329956 0.854901969432831]);
hold on;
plot(u42-U1d2,'LineWidth',3,'Color' ,[0 0 0]);
hold on;
plot(up2-U1d2,'LineWidth',3,'Color' ,[0.0784313753247261 0.168627455830574 0.549019634723663]);
set(gca,'linewidth',3);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20);

subplot(1,3,3)
plot(uf3-U1d3,'LineWidth',3,'Color',[0.388235300779343 0.545098066329956 0.854901969432831]);
hold on;
plot(u43-U1d3,'LineWidth',3,'Color' ,[0 0 0]);
hold on;
plot(up3-U1d3,'LineWidth',3,'Color' ,[0.0784313753247261 0.168627455830574 0.549019634723663]);
set(gca,'linewidth',3);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20);