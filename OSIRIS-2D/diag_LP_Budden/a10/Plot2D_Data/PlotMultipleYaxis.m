close all
% ------------------your 6 functions-----------------------------------------
x=0:0.5:10;
y1=0.1*sin(0.1*x);
y2=0.1*cos(0.1*x);
y3=0.1*cos(0.2*x);
y4=0.2*sin(0.4*x);
y5=0.7*sin(0.6*x);
y6=0.1*sin(0.8*x);
%-----------------------------------------------------------------------
xlab='time' ;              % x-axis title
pas=8;           % Number of ticks per axe
tail=8;             % font size
mark1='d';c1='g';    
mark2='<';c2='b';
mark3='s';c3='r';
mark4='p';c4='k';
mark5='o';c5='m';
mark6='*';c6='y';
fontname='courrier';
fontsize=9.2;
fontsize_ax=8;
% -----------your legend names--------------------------------------------
param1='tg\delta';
param2='Cx';
param3='Us';
param4='Q';
param5='L';
param6='H';
%-----------------y-axis titles----------------------------------------
paramm1='tg\delta .10^-^2';
paramm2='Cx(pF)';
paramm3='     Us(kV)';
paramm4='Q(pC)';
paramm5='L(mm)';
paramm6='H(m)';
%----------------------------------------------plot y1---------------------
plot(x,y1,c1)
ax1=gca;
set(gcf,'position',[100 100 1300 550])
line(x,y1,'color',c1,'Marker',mark1,'LineStyle','none','Markersize',tail,'parent',ax1);
set(ax1,'Ylim',double([min(y1) max(y1)]),'Xlim',[0 max(x)]);
ylabel(paramm1)
xlabel(xlab)
set(ax1,'position',[0.16 0.11 0.73 0.8],'fontsize',fontsize_ax,'Ycolor',c1)
pos=double(get(ax1,'position'));
%----------------------------------------------plot y2--------------------
ax2=axes('position',pos, 'XAxisLocation','bottom','YAxisLocation','left', 'Color','none',  'XColor',c2,'YColor',c2);
plot(x,y2,c2);
line(x,y2,'color',c2,'Marker',mark2,'LineStyle','none','Markersize',tail,'parent',ax2);
set(ax2,'Ylim',double([min(y2) max(y2)]),'Xlim',[0 max(x)],'Visible','off')
%----------------------------------------------plot y3-------------------
axe3=axes('position',pos, 'XAxisLocation','bottom','YAxisLocation','left','Color','none', 'XColor',c2,'YColor',c2);
plot(x,y3,c3);
line(x,y3,'color',c3,'Marker',mark3,'LineStyle','none','Markersize',tail,'parent',axe3);
set(axe3,'Ylim',double([min(y3) max(y3)]),'Xlim',[0 max(x)],'Visible','off')
%----------------------------------------------plot y4-----------------
axe4=axes('position',pos, 'XAxisLocation','bottom','YAxisLocation','left', 'Color','none', 'XColor',c2,'YColor',c2);
plot(x,y4,c4);
line(x,y4,'color',c4,'Marker',mark4,'LineStyle','none','Markersize',tail+2,'parent',axe4);
set(axe4,'Ylim',double([min(y4) max(y4)]),'Xlim',[0 max(x)],'Visible','off')
%----------------------------------------------plot y5--------------------
axe5=axes('position',pos, 'XAxisLocation','bottom','YAxisLocation','left', 'Color','none', 'XColor',c2,'YColor',c2);set(axe5,'Xlim',[0 max(x)])
plot(x,y5,c5);
line(x,y5,'color',c5,'Marker',mark5,'LineStyle','none','Markersize',tail,'parent',axe5);
set(axe5,'Ylim',double([min(y5) max(y5)]),'Xlim',[0 max(x)],'Visible','off');
%----------------------------------------------plot y6--------------------
axe6=axes('position',pos, 'XAxisLocation','bottom','YAxisLocation','left', 'Color','none', 'XColor',c2,'YColor',c2);set(axe5,'Xlim',[0 max(x)])
plot(x,y6,c6);
line(x,y6,'color',c6,'Marker',mark6,'LineStyle','none','Markersize',tail,'parent',axe6);
set(axe6,'Ylim',double([min(y6) max(y6)]),'Xlim',[0 max(x)],'Visible','off');
%---------------------------------------------------------ax12------------
pos12=[pos(1)+pos(3)-0.001 pos(2) 0.001 pos(4)]
axe12=axes('position',pos12,'XAxisLocation','top','YAxisLocation','right',...
'Color','none','fontsize',fontsize_ax, 'XColor',c2,'YColor',c2);
set(get(axe12,'XLabel'),'String',strvcat(' ',paramm2),'Fontname',fontname,'Fontsize',fontsize);
set(axe12,'Ylim',double([min(y2) max(y2)]))
inc2=abs(max(y2)-min(y2))/pas; 
my=min(y2):inc2:max(y2);
set(axe12,'Ytick',my)
%---------------------------------------------------------ax13-------------
pos13=[0.94 pos(2) 0.001 pos(4)]
axe13=axes('position',pos13, 'XAxisLocation','top','YAxisLocation','right',...
'Color','none','fontsize',fontsize_ax, 'XColor',c3,'YColor',c3);
set(get(axe13,'XLabel'),'String',strvcat(' ',paramm3),'Fontname',fontname,'Fontsize',fontsize)
set(axe13,'Ylim',double([min(y3)   max(y3)]))
inc3=(max(y3)-min(y3))/pas;
my=min(y3):inc3:max(y3);
set(axe13,'Ytick',my)
%---------------------------------------------------------ax14------------
pos14=[0.11 pos(2) 0.001 pos(4)]
axe14=axes('position',pos14, 'XAxisLocation','Bottom','YAxisLocation','left',...
'Color','none','fontsize',fontsize_ax,  'XColor',c4,'YColor',c4);
set(get(axe14,'XLabel'),'String',strvcat(' ',paramm4),'Fontname',fontname,'Fontsize',fontsize);
set(axe14,'Ylim',double([min(y4) max(y4)]))
inc4=(max(y4)-min(y4))/pas;
my=min(y4):inc4:max(y4);
set(axe14,'Ytick',my)
%---------------------------------------------------------ax15------------
pos15=[0.07 pos(2) 0.001 pos(4)]
axe15=axes('position',pos15,'XAxisLocation','bottom','YAxisLocation','left',...
'Color','none','fontsize',fontsize_ax, 'XColor',c5,'YColor',c5);
set(get(axe15,'XLabel'),'String',strvcat(' ',paramm5),'Fontname',fontname,'Fontsize',fontsize)
set(axe15,'Ylim',double([min(y5) max(y5)]))
inc5=(max(y5)-min(y5))/pas; 
my=min(y5):inc5:max(y5);
set(axe15,'Ytick',my)
%---------------------------------------------------------ax16------------
pos16=[0.03 pos(2) 0.001 pos(4)]
axe16=axes('position',pos16,'XAxisLocation','bottom','YAxisLocation','left',...
'Color','none','fontsize',fontsize_ax, 'XColor',c6,'YColor',c6);
set(get(axe16,'XLabel'),'String',strvcat(' ',paramm6),'Fontname',fontname,'Fontsize',fontsize)
set(axe16,'Ylim',double([min(y6) max(y6)]))
inc6=(max(y6)-min(y6))/pas; 
my=min(y6):inc6:max(y6);
set(axe16,'Ytick',my)
%------------------------------------------------legend------------------
ax20=axes('position',pos, 'Color','none')
line(-100,100,'color',c1,'Marker',mark1,'LineStyle','-','markerfacecolor',c1,'Markersize',tail,'parent',ax20);
hold on;line(-100,100,'color',c2,'Marker',mark2,'LineStyle','-','markerfacecolor',c2,'Markersize',tail,'parent',ax20);
hold on;line(-100,100,'color',c3,'Marker',mark3,'LineStyle','-','markerfacecolor',c3,'Markersize',tail,'parent',ax20);
hold on;line(-100,100,'color',c4,'Marker',mark4,'LineStyle','-','markerfacecolor',c4,'Markersize',tail+2,'parent',ax20);
hold on;line(-100,100,'color',c5,'Marker',mark5,'LineStyle','-','markerfacecolor',c5,'Markersize',tail,'parent',ax20);
hold on;line(-100,100,'color',c6,'Marker',mark6,'LineStyle','-','markerfacecolor',c6,'Markersize',tail,'parent',ax20);
set(ax20,'Xlim',[0 1]);
set(ax20,'visible','off');
grid(ax1);
name={param1;param2;param3;param4;param5;param6}
hleg=legend(ax20,name)
title(ax1,'figure1')

