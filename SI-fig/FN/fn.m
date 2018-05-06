%% SI-FN:The numerical solution to the FitzHugh-Nagumo equation

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file solve the FitzHugh-Nagumo equation presented in           %%%
%%% Ye Yuan, "Supplementary Information: Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%%% Method: The Finite Element Method                                   %%%
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used functions: 
%          pde_fn_ic: the initial condition
%          pde_fn_bc: the Boundary value condition
%          pde_fn_pde: the FitzHugh-Nagumo equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate mesh
m = 0;
x = linspace(-30,30,512); % space
t = linspace(0,100,401); % time

% Parameter 
A=0.2; B=0.002; C=0.001;
%A=0.2; B=0.0002; C=0.0001;
sol = pdepe(m,'pde_fn_pde','pde_fn_ic','pde_fn_bc',x,t,[],A,B,C);
u1 = sol(:,:,1);
u2 = sol(:,:,2);

% Plot
figure(1)
mesh(x,t,u1);
colormap('Winter')
xlabel('x')
ylabel('t')
zlabel('u(x,t)')
set(gca,'linewidth',2);
set(gca,'FontSize',14);
set(get(gca,'YLabel'),'Fontsize',14) 


figure(2)
mesh(x,t,u2);
colormap('Winter')
xlabel('x')
ylabel('t')
zlabel('w(x,t)')
set(gca,'linewidth',2);
set(gca,'FontSize',14);
set(get(gca,'YLabel'),'Fontsize',14) 
