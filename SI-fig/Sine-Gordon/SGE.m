%% SI-SG:The numerical solution to the sine-Gordon equation

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file solve the sine-Gordon equation presented in              %%%
%%% Ye Yuan, "Supplementary Information: Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%%% Method: the analytical solution                                     %%%
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used functions: No
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate mesh
x=linspace(-12.4,12.4,512);
t=linspace(0,20,256);
for i=1:256
for j=1:512
u1(i,j)=4*atan((sin(t(i)/sqrt(2)))/(cosh(x(j)/sqrt(2))));
end
end
mesh(t,x,u1')
colormap('Winter')
xlabel('t')
ylabel('x')
zlabel('u(x,t)')
%title('Numerical solution of Sine-Gordon equation')
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20) 

