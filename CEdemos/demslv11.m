% DEMSLV11 Generates figures for Chapter 3
function demslv11
disp('DEMSLV11 Generates figures for Chapter 3')
close all

% Function Iteration

xmin = 0.0;
xmax = 1.4;
xinit = 0.4;

g = inline('(x+0.2).^0.5');
xx = linspace(xmin,xmax);
yy = g(xx);

figure(1); 
plot(xx,yy,'LineWidth',2);
hold on
plot(xx,xx);
axis([min(xx) max(xx) min(xx) max(xx)])
title('Function Iteration')
axis square

x = xinit;
for it=1:20
   xnew = g(x);
   if it<4, plot([x x],[0 xnew],'k:'), end;  
   if it<4, plot([0 x],[xnew xnew],'k:'), end;
   plot([x x],[x xnew])  
   plot([x xnew],[xnew xnew])
   x = xnew;
end
plot([x x],[0.0 x],'k:')
fs=get(0,'DefaultTextFontSize')*(3/4);
h=text(x,0,'x^*'); 
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(0.40,0,'x^{(0)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(0.77,0,'x^{(1)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(0.99,0,'x^{(2)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(0.1,0,'45^o');
set(h,'VerticalAlignment','bottom','Fontsize',fs)
h=text(0,0.78,'x^{(1)}=g(x^{(0)})');
set(h,'HorizontalAlignment','right','Fontsize',fs)
h=text(0,0.99,'x^{(2)}=g(x^{(1)})');
set(h,'HorizontalAlignment','right','Fontsize',fs)
h=text(0,1.09,'x^{(3)}=g(x^{(2)})');
set(h,'HorizontalAlignment','right','Fontsize',fs)
set(gca,'XTick',[])
set(gca,'YTick',[])
hold off

% Newton's Method
 
xmin = 1.0;
xmax = 2.5;
xinit = 2.3;

xx = linspace(xmin,xmax);
yy = func(xx);
zz = zeros(size(xx));

figure(2); plot(xx,yy,'k',xx,zz,'k','LineWidth',2);
title('Newton''s Method')
axis square

hold on
x = xinit;
for it=1:4
   [y d] = func(x);
   xnew = x - y/d;
   plot([x x],[0 y],'k:')
   plot([x xnew],[y 0],'k')
   x = xnew;
end
h=text(2.30,0,'x^{(0)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(1.77,0,'x^{(1)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(1.42,0,'x^{(2)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(1.19,0,'x^*');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
set(gca,'XTick',[])
set(gca,'YTick',[0])
hold off


% Secant Method

figure(3); plot(xx,yy,'k',xx,zz,'k','LineWidth',2);
title('Secant Method')
axis square

hold on
x0 = xinit + 0.1;
x1 = xinit - 0.1;
y0 = func(x0);
y1 = func(x1);
plot([x0 x0],[0 y0],'k:')
plot([x1 x1],[0 y1],'k:')
h=text(x0,0,'x^{(0)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(x1,0,'x^{(1)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
for it=1:2
   x = x1-y1*(x1-x0)/(y1-y0);
   y = func(x);
   plot([x0 x],[y0 0],'k')
   plot([x x1],[0 y1],'k');
   plot([x x],[0 y],'k:')
   x0 = x1;
   y0 = y1;
   x1 = x;
   y1 = y;
end
h=text(1.76,0,'x^{(2)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(1.52,0,'x^{(3)}');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
h=text(1.19,0,'x^*');
set(h,'VerticalAlignment','top','HorizontalAlignment','center','Fontsize',fs)
set(gca,'XTick',[])
set(gca,'YTick',[0])
hold off



prtfigs(mfilename,'Function Iteration',1)

prtfigs(mfilename,'Newton''s Method',2)

prtfigs(mfilename,'Secant Iteration',3)


%%%%%%%%%%%%%%%%%

function [y,d] = func(x);
alpha = 4;
y = x.^alpha-2;
d = alpha*x.^(alpha-1);
