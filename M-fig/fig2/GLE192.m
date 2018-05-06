%%M_Fig_2.The Schematics of the S^3d method

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file generate the Fig2-M presented in                          %%%
%%% Ye Yuan, " Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;

%% download cal/car06212
load('data192.mat')
R = real(AR);
I = imag(AR);
AR = R-1i*I;


%% Plot 3-D and contour
figure(1)
B=abs(AR)+abs(AL);
B=B';
figure(1)
C=B+0.02;
figure(1)
mesh(C(1:200,:))
colormap('Winter')
hold on
contour(C(1:200,:))
axis off

%% make grid
% space
L=82.0119; 
n=180;
x2=-L/2:0.458167:L/2;
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
% time
t=0:1.5625:1187; 


%% download the identified coefficients of S^3d
load('c192.mat');

% physical parameters
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
I=658; 
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




%% reconstruction 192 

figure(2)
subplot(1,2,1),mesh(abs(AR(:,1:760))+abs(AL(:,1:760)));
colormap('Winter')
axis([0 760 0 180]);
set(gca,'XTick',[0,200,400,600,760], 'YTick',[0, 30, 60, 90,120,150, 180]);
xlabel('t')
ylabel('x')
%zlabel('A_R(x,t)+A_L(x,t)')
%title('Experiment','FontSize',40)
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',30) 
set(get(gca,'XLabel'),'Fontsize',30) 
% angle={[0 90]};
% view(angle{1})
subplot(1,2,2),mesh(abs(usol(:,1:760))+abs(vsol(:,1:760)));
colormap('Winter')
axis([0 760 0 180]);
set(gca,'XTick',[0,200,400,600,760], 'YTick',[0, 30, 60, 90,120,150, 180]);
xlabel('t')
ylabel('x')
%zlabel('A_R(x,t)+A_L(x,t)')
%title('Reconstruction','FontSize',40)
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',30) 
set(get(gca,'XLabel'),'Fontsize',30) 
% angle={[0 90]};
% view(angle{1})


%% plot initial value
figure(3)
plot(abs(AL(:,I)),'LineWidth',4,'Color',[0 0.200000002980232 0.600000023841858]);
   axis([0 180 0 0.02]);
set(gca,'XTick',[0, 30, 60, 90,120,150, 180], 'YTick',[0,0.005,0.01,0.015,0.02]);
    hold on
 plot(abs(AR(:,I)),'LineWidth',4,'Color',[0 0.600000023841858 0]);
   axis([0 180 0 0.02]);
set(gca,'XTick',[0, 30, 60, 90,120,150, 180], 'YTick',[0,0.005,0.01,0.015,0.02]);
set(gca,'linewidth',2);
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Fontsize',30) 
annotation(figure(3),'textbox',...
    [0.52783143963729 0.775375239969119 0.241038216560509 0.121272365805162],...
    'String',{'cal/car06192'},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');


%% plot deritive A_t 
figure(4)
load('A_dt.mat')
plot(abs(AL_dt(1,21:780)),'LineWidth',3,...
    'Color',[0.635294139385223 0.0784313753247261 0.184313729405403])

