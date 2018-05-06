%%SI_car/cal06182. Reconstruction 

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file generate the Fig.15-SI presented in                                %%%
%%% Ye Yuan, " Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;

%% download cal/car06182
load('data182.mat')
R = real(AR);
I = imag(AR);
AR = R-1i*I;
% space
L=82.0119; 
n=180;
x2=-L/2:0.458167:L/2;
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
% time
t=0:1.5625:312.5;


%% download the identified coefficients of S^3d
load('c182.mat');

ralp=wr(1);
repstau=wr(2);
rxi2tau=wr(3);
rgtau=wr(4);
rhtau=wr(5);
repsc0tau=wr(6);
rc1xi2tau=wr(7);
rgc2tau=wr(8);
rhc3tau=wr(9);
rbeta0=(rxi2tau+rc1xi2tau*1i);
rbeta1=(repstau+repsc0tau*1i);
rbeta2=(rgtau+rgc2tau*1i);
rbeta3=(rhtau+rhc3tau*1i);

lalp=wl(1);
lepstau=wl(2);
lxi2tau=wl(3);
lgtau=wl(4);
lhtau=wl(5);
lepsc0tau=wl(6);
lc1xi2tau=wl(7);
lgc2tau=wl(8);
lhc3tau=wl(9);
lbeta0=(lxi2tau+lc1xi2tau*1i);
lbeta1=(lepstau+lepsc0tau*1i);
lbeta2=(lgtau+lgc2tau*1i);
lbeta3=(lhtau+lhc3tau*1i);

% initial conditions
I=1;
u=AR(:,I);
v=AL(:,I);

% solve 
ut=fft(u);
vt=fft(v);
uvt=[ut vt];
[t,uvsol]=ode45('gle_rhs',t,uvt,[],rbeta0,rbeta1,rbeta2,rbeta3,ralp,lbeta0,lbeta1,lbeta2,lbeta3,lalp,k,n);
for j=1:length(t)
usol(j,:)=ifft(uvsol(j,1:n)); 
vsol(j,:)=ifft(uvsol(j,n+1:2*n));
end
usol=usol.';
vsol=vsol.';

%% Reconstruction the dynamics underlying data (car/cal06182)
figure(1)
subplot(1,2,1),mesh(abs(AR(:,1:200))+abs(AL(:,1:200)));
colormap('Winter')
axis([0 200 0 180]);
set(gca,'XTick',[0,50,100,150,200], 'YTick',[0, 30, 60, 90,120,150, 180]);
%xlabel('t')
ylabel('x')
zlabel('A_R(x,t)+A_L(x,t)')
%title('Experiment','FontSize',40)
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',30) 
set(get(gca,'XLabel'),'Fontsize',30) 
angle={[0 90]};
view(angle{1})
subplot(1,2,2),mesh(abs(usol(:,1:200))+abs(vsol(:,1:200)));
colormap('Winter')
axis([0 200 0 180]);
set(gca,'XTick',[0,50,100,150,200], 'YTick',[0, 30, 60, 90,120,150, 180]);
%xlabel('t')
ylabel('x')
zlabel('A_R(x,t)+A_L(x,t)')
%title('Reconstruction','FontSize',40)
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',30) 
set(get(gca,'XLabel'),'Fontsize',30) 
angle={[0 90]};
view(angle{1})

%% plot initial value
figure(2)
    plot(abs(AL(:,I)),'LineWidth',4,'Color',[0 0.200000002980232 0.600000023841858]);
   axis([0 0.02 0 180]);
set(gca,'XTick',[0,0.005,0.01,0.015,0.02], 'YTick',[0, 30, 60, 90,120,150, 180]);
    hold on
   plot(abs(AR(:,I)),'LineWidth',4,'Color',[0 0.600000023841858 0]);
   axis([0 180 0 0.02]);
set(gca,'XTick',[0, 30, 60, 90,120,150, 180], 'YTick',[0,0.005,0.01,0.015,0.02]);
xlabel('space')
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',30) 
annotation(figure(2),'textbox',...
    [0.52783143963729 0.775375239969119 0.241038216560509 0.121272365805162],...
    'String',{'cal/car06182'},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% plot
% figure(3)
% subplot(2,2,1),mesh(abs(usol));
% colormap('Winter')
% subplot(2,2,3),mesh(abs(vsol));
% colormap('Winter')
% subplot(2,2,2),mesh(abs(AR(:,1:400)));
% colormap('Winter')
% subplot(2,2,4),mesh(abs(AL(:,1:400)));
% colormap('Winter')