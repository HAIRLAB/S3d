%% SI-KS: The numerical solution to the KS equation

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file plot the dynamics of the KS equation presented in         %%%
%%% Ye Yuan, "Supplementary Information: Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used functions: No
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('KS.mat');
u=uu;


% Plot the results
figure(1)
mesh(tt,x, u); % 3-D surface plot
colormap('Winter')
xlabel('t')
ylabel('x')
zlabel('u(x,t)')
set(gca,'linewidth',2);
set(gca,'FontSize',14);
set(get(gca,'YLabel'),'Fontsize',14) 