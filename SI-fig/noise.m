%% FigS_FLOW_GEN. A flow diagram for identfying the PDE

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file generate the FigS1 presented in                           %%%
%%% Ye Yuan, "Supplementary Information: Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clean data
load('kdv.mat');

%% Noise data
eps = 0.01;
u_std = std(reshape(u,1,size(u,1)*size(u,2)));
u_e=a + eps*u_std*randn(size(u));


%% Show results
figure(1),clf,hold on;
plot(x,u(:,1),'MarkerSize',0.5,'LineWidth',3,'Color',[0.800000011920929 0 0]);
hold on
plot(x,u_e(:,1),'MarkerSize',15,'Marker','.','LineWidth',2.5,'LineStyle','none',...
    'Color',[0.0784313753247261 0.168627455830574 0.549019634723663]);
% plot(x,u_e(:,500),'MarkerSize',0.5,'LineWidth',2.5,'LineStyle','-.',...
%     'Color',[0.800000011920929 0 0])
% plot(x,u_e(:,1000),'MarkerSize',0.5,'LineWidth',2.5,'LineStyle','-.',...
%     'Color',[0 0 0])
%legend('t=1','t=500','x=1000','Location','Southeast')
xlabel('x')
ylabel('snapshot')
title([ 'the sanpshot randomly sampled from the grid']);
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20) 



%%
figure(2)
plot(x,u(:,1)-u_e(:,1),'MarkerSize',0.5,'LineWidth',2.5,'LineStyle','-.',...
    'Color',[0.800000011920929 0 0]);
xlabel('t')
ylabel('error')
title([ 'the error with x=200']);
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20) 
%%
figure(3)
plot(x,u(:,500)-u_e(:,500),'MarkerSize',0.5,'LineWidth',2.5,'LineStyle','-.',...
    'Color',[0 0 0])
xlabel('t')
ylabel('error')
title([ 'the error with x=212']);
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20) 


%% Subsample 
[m,n]=size(u);
for i=1:m
er(i,1)=norm(u(:,i),2);
end
err=1/sqrt(512*501)*norm(er(:,1),2);


%% Subsample 
figure(4)
for i=51:1:100
plot(t,u(i,:),'MarkerSize',0.5,'LineWidth',2.5,'LineStyle','-');hold on
end
xlabel('time')
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20) 