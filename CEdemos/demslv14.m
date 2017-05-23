% DEMSLV14 Illustrates Problematic Complementarity Problems
disp('DEMSLV14 Illustrates Problematic Complementarity Problems')
close all

xmin = -1.5;
xmax =  1.5;
ymin = -5;
ymax =  5;

x = nodeunif(500,xmin,xmax);
y = -4*x;
z = zeros(size(x));

close all


fs=get(0,'DefaultTextFontSize')*(3/4);

vtext=-0.2;

figure(1)

subplot(2,2,1)
hold on
title('f''<0, f(a)<0')
a = 1;
b = 3;
ym = minmax(x,a,b,y);
hold on
plot(x,ym,'LineWidth',4);
plot(x,y ,'LineWidth',2);
plot(x,z , ':');
plot(x,a-x,'--');
plot(x,b-x,'--');

plot(a,0,'*');
plot(b,0,'*');
h=text(a,vtext,'a');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
set(gca,'Xtick',[ ])
set(gca,'Ytick',[0])
axis square
axis ([xmin xmax ymin ymax])
hold off

subplot(2,2,2)
hold on
title('f''<0, f(b)>0')
a = -3;
b = -1;
ym = minmax(x,a,b,y);
hold on
plot(x,ym,'LineWidth',4);
plot(x,y ,'LineWidth',2);
plot(x,z , ':');
plot(x,a-x,'--');
plot(x,b-x,'--');
plot(b,0,'*')
h=text(b,vtext,{'b'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
set(gca,'Xtick',[ ])
set(gca,'Ytick',[0])
axis square
axis ([xmin xmax ymin ymax])
hold off

subplot(2,2,3)
hold on
title('f''<0, f(a)>0>f(b)')
a = -1;
b =  1;
ym = minmax(x,a,b,y);
hold on
plot(x,ym,'LineWidth',4);
plot(x,y ,'LineWidth',2);
plot(x,z , ':');
plot(x,a-x,'--');
plot(x,b-x,'--');
plot(a,0,'*');plot( 1,0,'*')
h=text(a,vtext,{'a'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(b,vtext,{'b'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
set(gca,'Xtick',[ ])
set(gca,'Ytick',[0])
axis square
axis ([xmin xmax ymin ymax])
hold off

subplot(2,2,4)
hold on
title('f''>0')
y  = -y;
a = -1;
b =  1;
ym = minmax(x,a,b,y);
hold on
plot(x,ym,'LineWidth',4);
plot(x,y ,'LineWidth',2);
plot(x,z , ':');
plot(x,a-x,'--');
plot(x,b-x,'--');
plot(-1,0,'*')
plot( 1,0,'*')
h=text(a,vtext,{'a'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(b,vtext,{'b'});
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
set(gca,'Xtick',[ ])
set(gca,'Ytick',[0])
axis square
axis ([xmin xmax ymin ymax])
hold off

prtfigs(mfilename,'Reformulated Complementarity Problems',1)
