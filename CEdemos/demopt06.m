% DEMOPT06 Illustrates constrained optimization problems
function demopt06
close all

disp(' ')
disp('DEMOPT06 Illustrates constrained optimization problems')

x=nodeunif(100,-0.5,1.5);
a=0;
b=1;

close all
figure(1)
f=inline('1.5-(x-.75).^2','x');
subplot(1,2,1)
plot(x,f(x));
hold on;
plot([a;a],[-0.5;2],'--');
plot([b;b],[-0.5;2],'--');
axis([-0.5 1.5 -0.5 2])
axis square
title('f(x) = 1.5-(x-3/4)^2, x^* = 3/4')
hold off


f=inline('(x-.75).^2','x');
subplot(1,2,2)
plot(x,f(x));
hold on;
plot([a;a],[-0.5;2],'--');
plot([b;b],[-0.5;2],'--');
axis([-0.5 1.5 -0.5 2])
axis square
title('f(x) = -2+(x-3/4)^2, x^* = 0 & 1')
hold off


figure(2)
f=inline('-2*(x-3/4)','x');
subplot(1,2,1)
y  = f(x);
ym = minmax(x,a,b,y);
plot(x,a-x,'--');
hold on;
plot(x,b-x,'--');
plot(x,y);
plot(x,ym,'LineWidth',2);
plot(x,zeros(length(x),1),':');
axis([-0.5 1.5 -2 2])
axis square
title('f''(x) = -2(x-3/4)')
hold off


f=inline('2*(x-0.75)','x');
subplot(1,2,2)
y  = f(x);
ym = minmax(x,a,b,y);
plot(x,a-x,'--');
hold on;
plot(x,b-x,'--');
plot(x,y);
plot(x,ym,'LineWidth',2);
plot(x,zeros(length(x),1),':');
axis([-0.5 1.5 -2 2])
axis square
title('f''(x) =  2(x-3/4)')
hold off

prtfigs(mfilename,'One-Dimensional Maximization Problems',1)
prtfigs(mfilename,'Complementarity Conditions for Maximization Problems',2)
