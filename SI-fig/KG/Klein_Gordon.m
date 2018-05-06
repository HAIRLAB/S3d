%% SI-KG:The numerical solution to the Klein-Gordon equation

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file solve the Klein-Gordon equation presented in              %%%
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

xl=0; xr=1;       % x domain [xl,xr]
J = 100;          % J: number of division for x
dx = (xr-xl) / J; % dx: mesh size
tf = 3;           % final simulation time
Nt = 1000;        % Nt: number of time steps
dt = tf/Nt;
mu = dt/(dx)^2;

% Parameter
a = -1;
b = 1;
r = 1;
B = b/r;
c = 0.5;
K = b/r;

% Evaluate the initial conditions
x = xl : dx : xr; % generate the grid point
% f(1:J+1) since array index starts from 1


% store the solution at all grid points for all time steps
u = zeros(J+1,Nt+1);

% Find the approximate solution at each time step
for n = 1:Nt
    % boundary condition at left side
    gl = 0;
    % boundary condition at right side
    gr = 0;
    if n==1 % first time step
        u(:,1) = B*K*x.^2;      
    elseif n==2
        for j=2:J % interior nodes
            u(j,2) = u(j,1) + dt*B*c*K*2*2*dt;
        end
        u(1,n) = B*K*x(1);            % the left-end point
        u(J+1,n) = B*K*x(J+1);        % the right-end point      
    else      
        for j=2:J                    % interior nodes
            u(j,n)=2*u(j,n-1)-u(j,n-2)+dt^2*(-a*( u(j+1,n-1)-2*u(j,n-1)+...
                u(j-1,n-1))/dx^2 - b*u(j,n-1) - r*u(j,n-1)^3);
        end
        u(1,n) = B*K*x(1)^2;         % the left-end point
        u(J+1,n) = B*K*x(J+1)^2;     % the right-end point
        
    end


end
% Plot the results
tt = dt:dt:(Nt+1)*dt;
figure(1)
mesh(x,tt, u'); % 3-D surface plot
colormap('Winter')
xlabel('x')
ylabel('t')
zlabel('u(x,t)')
title('Numerical solution of Klein Gordon equation')
set(gca,'linewidth',2);
set(gca,'FontSize',14);
set(get(gca,'YLabel'),'Fontsize',14) 
