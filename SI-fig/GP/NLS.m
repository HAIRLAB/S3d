%% SI-GP:The numerical solution to the GP equation

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file solve the GP equation presented in                        %%%
%%% Ye Yuan, "Supplementary Information: Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%%% Method: The spectral Method                                         %%%
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used functions: 
%          nls_rhs: equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate mesh

L=2*pi; n=512;
x2=linspace(-L/2,L/2,n+1); x=x2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
% time
t=linspace(0,8,501);
% initial conditions
N=2;
u = 4*x.^2.*exp(-2.*x.^2).*exp(1i);
ut=fft(u);
[t,utsol]=ode45('nls_rhs',t,ut,[],k);
for j=1:length(t)
usol(j,:)=ifft(utsol(j,:)); % bring back to space
end
usol=usol.';

%% plot
mesh(x,t,abs(usol));
colormap('Winter')
xlabel('x')
ylabel('t')
zlabel('u(x,t)')
set(gca,'linewidth',2);
set(gca,'FontSize',14);
set(get(gca,'YLabel'),'Fontsize',14) 