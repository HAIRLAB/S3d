%%SI_car/cal06252. Structural analysis

%%%%%%%%%%%%%%%%%%%%%%%%%  S^3d method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file generate the Fig.14-SI presented in                       %%%
%%% Ye Yuan, " Machine Discovery of Partial Differential 
%%%             Equations from Spatiotemporal Data"                     %%%
%   @author: Ye Yuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('fig252left.mat');
x = num;
y1 = number;
y2 = err;
subplot(1,2,1);
[AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
set(AX(1),'XColor','k','YColor',[0 0 0]);
set(AX(2),'XColor','k','YColor','w');
set(gca,'color',[0.941176474094391 0.941176474094391 0.941176474094391]);
set(gcf,'color','w');
HH1=get(AX(1),'Ylabel');
set(HH1,'String','Number of terms');
set(HH1,'color',[0 0 0]);
HH2=get(AX(2),'Ylabel');
set(HH2,'color','r');
set(H1,'MarKer','o','Linestyle','none','MarkerFaceColor',[0 0 0.600000023841858],...
    'Color',[0 0 0.600000023841858],'MarkerSize',4); 
set(AX(1),'ylim',[0,20],'ytick',0:5:20)%������y�����ݷ�Χ��0��0.5�����̶���ʾ��0,0.1,0.2...0.5��  
set(H2,'Linestyle','-','color',[0.635294139385223 0.0784313753247261 0.184313729405403],'Linewidth',3);
set(AX(2),'ylim',[0.16,0.32],'ytick',0.16:0.02:0.32)%������y�����ݷ�Χ��0��3�����̶���ʾ��0,1,2,3��  
xlabel('\lambda');
title('cal06252');
set(gca,'linewidth',2);
set(gca,'FontSize',25);
set(get(gca,'YLabel'),'Fontsize',35);

load('fig252right.mat');
x = num;
y1 = number;
y2 = err;
subplot(1,2,2);
[AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
set(AX(1),'XColor','k','YColor',[0 0 0]);
set(AX(2),'XColor','k','YColor','w');
set(gca,'color',[0.941176474094391 0.941176474094391 0.941176474094391]);
set(gcf,'color','w');
HH1=get(AX(1),'Ylabel');
set(HH1,'String','Number of terms');
set(HH1,'color',[0 0 0]);
HH2=get(AX(2),'Ylabel');
set(HH2,'color','r');
set(H1,'MarKer','o','Linestyle','none','MarkerFaceColor',[0 0 0.600000023841858],...
    'Color',[0 0 0.600000023841858],'MarkerSize',4); 
set(AX(1),'ylim',[0,20],'ytick',0:5:20)%������y�����ݷ�Χ��0��0.5�����̶���ʾ��0,0.1,0.2...0.5��  
set(H2,'Linestyle','-','color',[0.635294139385223 0.0784313753247261 0.184313729405403],'Linewidth',3);
set(AX(2),'ylim',[0.20,0.38],'ytick',0.20:0.02:0.38)%������y�����ݷ�Χ��0��3�����̶���ʾ��0,1,2,3��  
xlabel('\lambda');
title('car06252');
set(gca,'linewidth',2);
set(gca,'FontSize',25);
set(get(gca,'YLabel'),'Fontsize',35);