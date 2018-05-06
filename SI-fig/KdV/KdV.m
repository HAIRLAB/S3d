%-----------------------------------------------------
% soliton.m: solve the double solition problem
%
% Functions used:
% ushu51.m: the exact solution
% resshu5.m: define the residual function
% reconuxxxp.m: reconstruct u_{xxx} from u
% reconuxp.m: reconstruct u_x from u
% filterLI.m: our 10th order filter
%-----------------------------------------------------
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate 1D mesh
% for ex 4.5 (4.9) & (4.10)
Nx=256; XL=0; XR=2;
% for ex 4.5 (4.11)
%Nx=151; XL=0; XR=3;
dx=(XR-XL)/(Nx-1);
for i=1:Nx
XX(i)=XL+(i-1)*dx;
end
%%%%%%%%%march-in-time%%%%%%%%%%%%
% cfl has to be small!
cfl=0.05; dt=cfl*dx^2; tlast=2;
Nstep=tlast/dt;
IRK=2;
%ap=5*10^(-4); % coefficient eps for the PDE, (4.9)
ap=4.84*10^(-4); % for (4.10) of Shu
%ap=10^(-4); % for (4.11)
% calculate init cond.: stored at u0ex
[u0ex]=feval(@ushu51,XX,ap);
% for ex.4.1., ap likes t
%[u0ex]=feval(@ushu51,XX,0);
uold=u0ex;
IT=1; % indicator for storing solution uall(1:Nx,1:NT)
uall(:,1)=u0ex(:); % for contour plot
% we use 4th-order RK scheme:IRK=2,4
for k=1:Nstep
if IRK==4
% u0x=reconux(uold,Nx,dx); % reconstruct ux from u_n
% k0=dt*resshu5(u0x); % k0=dt*R((u_n)_x)
% u1=uold+k0/2;
% u1x=reconux(u1,Nx,dx);
% k1=dt*resshu5(u1x);
% u2=u1+k1/2;
% u2x=reconux(u2,Nx,dx);
% k2=dt*resshu5(u2x);
% u3=u2+k2;
% u3x=reconux(u3,Nx,dx);
% k3=dt*resshu5(u3x);
% unew=uold+(k0+2*k1+2*k2+k3)/6.; % finish one-time step
elseif IRK==2
% method 1: use 1st derivative construction three times
%u0x=reconuxp(uold,Nx,dx); % reconstruct ux from u_n
%u0xx=reconuxp(u0x,Nx,dx);
%u0xxx=reconuxp(u0xx,Nx,dx); % obtain u_xxx
% method 2: construct directly
u0xxx=reconuxxxp(uold,Nx,dx);
u2x=reconuxp(uold.^2,Nx,dx); % reconstruct u^2
k0=dt*resshu5(u2x,u0xxx,ap); % k0=dt*R((u_n)_x)
u1=uold+k0;
%u1x=reconuxp(u1,Nx,dx);
%uxx=reconuxp(u1x,Nx,dx);
%uxxx=reconuxp(uxx,Nx,dx); % obtain u_xxx
uxxx=reconuxxxp(u1,Nx,dx);
u2x=reconuxp(u1.^2,Nx,dx); % reconstruct u^2
k1=dt*resshu5(u2x,uxxx,ap);
u2=u1+k1;
unew=(uold+u2)/2;
end
uold=unew; % update for next time step
% for (4.11)
subplot(2,2,1), plot(XX,u0ex,'g-');
axis([XL XR -0.1 1]); xlabel('x'); ylabel('u(x,t=0)');
if abs(k*dt-1) < eps
    subplot(2,2,2),
unew=filterLI(unew,Nx);
plot(XX,unew,'g-'); % plot t=1
axis([XL XR -0.1 1]); xlabel('x'); ylabel('u(x,t=1)');
elseif abs(k*dt-2) < eps
subplot(2,2,3),
unew=filterLI(unew,Nx);
plot(XX,unew,'g-'); % plot t=2
axis([XL XR -0.1 1]); xlabel('x'); ylabel('u(x,t=2)');
end
% save some vaules for contour plot
if mod(k,500)==0
disp(k), IT=IT+1;
uall(:,IT)=unew(:);
end
end
%% do contour plot
ht=tlast/(IT-1);
[X,Y]=meshgrid(XL:dx:XR, 0:ht:tlast);
mesh(Y,X,uall');
axis([XL XR 0 tlast -0.1 1]);
xlabel('t')
ylabel('x')
zlabel('u(x,t)')
%title('Numerical solution of Sine-Gordon equation')
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',20) 


