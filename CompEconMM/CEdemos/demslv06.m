%% DEMSLV06 Nonlinear Equation Methods
% 
%  Illusttates function iteration, Newton, and secant methods graphically.

function demslv06

% Preliminary tasks
demosetup(mfilename)


%% Function Iteration

xmin = 0.0;
xmax = 1.4;
xinit = 0.3;
xstar = 0.5*(1+sqrt(1.8));
xx = linspace(xmin,xmax);
yy = g(xx);
x = zeros(20,1);
y = zeros(20,1);

figure
hold on
plot(xx,yy)
plot(xx,xx,'k-','LineWidth',2)
axis([min(xx) max(xx) min(xx) max(xx)])
title('Function Iteration')
x(1) = xinit;
y(1) = g(x(1));
for i=2:20
  x(i) = y(i-1);
  y(i) = g(x(i));
end
plottext(x(1),0,'$x_0$','center','top')
plottext(x(2),0,'$x_1$','center','top')
plottext(x(3),0,'$x_2$','center','top')
plottext(xstar,0,'$x*$','center','top')
plottext(0,x(2),'$x_1$','right','middle')
plottext(0,x(3),'$x_2$','right','middle')
plottext(0,x(4),'$x_3$','right','middle')
plottext(0,xstar,'$x^*$','right','middle')
for i=[1:3 19]
  plot([x(i) x(i)],[0 x(i+1)],'k:','LineWidth',2);
  plot([0 x(i)],[x(i+1) x(i+1)],'k:','LineWidth',2);
end
for i=2:20
  plot([x(i-1) x(i-1)],[x(i-1) x(i)],'r','LineWidth',2)
  plot([x(i-1) x(i)],[x(i) x(i)],'r','LineWidth',2)
end
for i=1:3
  plotbullet(x(i),x(i))
  plotbullet(x(i),x(i+1))
end
plottext(0.15,0.01,'$45^\circ$','center','bottom',12)
plottext(xmax,sqrt(xmax+.2),'$g$','center','top')
plottext(xmin-.02,0,'0','right','middle')
plotbullet(xstar,xstar)
set(gca,'XTick',[])
set(gca,'YTick',[])


%% Newton's Method

xmin = 1.0;
xmax = 2.5;
xinit = xmax-0.1;
xstar = 3^(1/5);
xx = linspace(xmin,xmax);
yy = f(xx);
zz = zeros(size(xx));

figure
hold on
plot(xx,yy)
plot(xx,zz,'k')
title('Newton''s Method')
axis off

x(1) = xinit;
y(1) = f(x(1));
for i=2:4
  [ylag,dlag] = f(x(i-1));
  x(i) = x(i-1) - ylag/dlag;
  y(i) = f(x(i));
end
plottext(x(1),0,'$x_0$','center','top')
plottext(x(2),0,'$x_1$','center','top')
plottext(x(3),0,'$x_2$','center','top')
plottext(x(4),0,'$x_3$','center','top')
plottext(xstar,0,'$x^*$','center','top')
for i=2:4
  plot([x(i-1) x(i)],[y(i-1) 0],'r','LineWidth',2)
end
for i=1:4
  plot([x(i) x(i)],[0 y(i)],'k:','LineWidth',2)
  plotbullet(x(i),y(i))
  plotbullet(x(i),0)
end
plotbullet(xstar,0)
plottext(xmin-.02,0,'0','right','middle')
plottext(xmax,f(xmax),'$f$','center','bottom')
set(gca,'XTick',[])
set(gca,'YTick',[])


%% Secant Method

figure
hold on
plot(xx,yy);
plot(xx,zz,'k');
title('Secant Method')
axis off

x(1) = xinit;
x(2) = xinit-0.25;
y(1) = f(x(1));
y(2) = f(x(2));
for i=3:4
  x(i) = x(i-1)-y(i-1)*(x(i-1)-x(i-2))/(y(i-1)-y(i-2));
  y(i) = f(x(i));
end
plottext(x(1),0,'$x_0$','center','top')
plottext(x(2),0,'$x_1$','center','top')
plottext(x(3),0,'$x_2$','center','top')
plottext(x(4),0,'$x_3$','center','top')
plottext(xstar,0,'$x^*$','center','top')
for i=3:4
  plot([x(i) x(i-2)],[0 y(i-2)],'r','LineWidth',2);
end
for i=1:4
  plot([x(i) x(i)],[0 y(i)],'k:','LineWidth',2)
  plotbullet(x(i),y(i))
  plotbullet(x(i),0)
end
plotbullet(xstar,0)
plottext(xmin-.02,0,'0','right','middle')
plottext(xmax,f(xmax),'$f$','center','bottom')
set(gca,'XTick',[])
set(gca,'YTick',[])


%% SAVE FIGURES
printfigures(mfilename)


%% Function
function [fval,fjac] = f(x)
fval = x.^5-3;
fjac = 5*x.^4;

%% Function
function gval = g(x)
gval = (x+0.2).^0.5';