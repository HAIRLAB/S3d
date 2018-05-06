%% SI-KG:The numerical solution to the Fisher equation

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file solve the Fisher equation presented in                    %%%
%%% Ye Yuan, "Supplementary Information: Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%%% Method: The Finite Difference Method                                %%%
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used functions: No
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate mesh
xl=-6; xr=6;         % x domain [xl,xr]
J = 200;             % J: number of division for x
dx = (xr-xl) / J;    % dx: mesh size
tf = 10;             % final simulation time
Nt = 1000;           % Nt: number of time steps
dt = tf/Nt;
mu = dt/(dx)^2;
t=linspace(0,10,1000);

% Evaluate the initial conditions
x = xl : dx : xr; % generate the grid point
% f(1:J+1) since array index starts from 1
for i = 1:J+1
    xx = ((i-1)*dx+xl);
    if xx<-1
        f(i) = exp(10*(xx+1));
    elseif xx>1
        f(i) = exp(-10*(xx-1));
    else f(i)=1;
        
    end
end
v = 0.1;

% store the solution at all grid points for all time steps
u = zeros(J+1,Nt);

% Find the approximate solution at each time step
for n = 1:Nt
    
    % boundary condition at left side
    gl = 0;
    % boundary condition at right side
    gr = 0;
    if n==1 % first time step
        for j=2:J % interior nodes 
            u(j,1)=(0.1*(f(j+1)-2*f(j)+f(j-1))/dx^2 + f(j)-f(j)^2)*dt+f(j);
        end
        u(1,n) = gl; % the left-end point
        u(J+1,n) = gr; % the right-end point

    else
        for j=2:J % interior nodes
             u(j,n) =( 0.1*( u(j+1,n-1) - 2*u(j,n-1) +u(j-1,n-1) )/dx^2+...
                 u(j,n-1) - u(j,n-1)^2 )*dt + u(j,n-1);
        end
        u(1,n) = gl; % the left-end point
        u(J+1,n) = gr; % the right-end point

    end
    % calculate the analytic solution
    for j=1:J+1
        u_exa(j,n)=(1/(1+exp(sqrt(1/6))*x(j)-(5/6)*t(n))).^2;
    end
    
end
%% Plot 3_D mesh 
figure(1)
mesh(x,t, u'); % 3-D surface plot
colormap('Winter')
xlabel('x')
ylabel('t')
zlabel('u(x,t)')
title('Numerical solution of Fisher equation')
set(gca,'linewidth',2);
set(gca,'FontSize',14);
set(get(gca,'YLabel'),'Fontsize',14) 

%% Plot the  countor
figure(2)
A=u';
for j=100:30:300
plot(x,A(j,:),'MarkerSize',0.5,'LineWidth',2.5,'LineStyle','-.');
hold on
end
hold on
for j=1:6:40
plot(x,A(j,:),'MarkerSize',0.5,'LineWidth',2.5,'LineStyle','-.');
hold on
end
hold on
for j=500:80:1000
plot(x,A(j,:),'MarkerSize',0.5,'LineWidth',2.5,'LineStyle','-.');
hold on
end
xlabel('x')
ylabel('u(x,t)')
set(gca,'linewidth',2);
set(gca,'FontSize',14);
set(get(gca,'YLabel'),'Fontsize',14) 
