%% DEMMATH05 Inverse Function, Implicit Function and Mean Value Theorems

% Preliminary tasks
demosetup(mfilename)
optset('broyden','showiters',0)


%% Inverse Function Theorem

f  = @(p)  0.50*p.^-0.2 + 0.500*p.^-0.5 + 1;
f1 = @(p) -0.10*p.^-1.2 - 0.250*p.^-1.5;
f2 = @(p)  0.12*p.^-2.2 + 0.375*p.^-2.5;

n  = 200;
p  = nodeunif(n,0.4,3.0);
q  = f(p);
g1 = 1/f1(1);
g2 = -f2(1)/f1(1)^3;
p1 = 1+g1*(q-2);
p2 = 1+g1*(q-2) + 0.5*g2*(q-2).^2;

figure
plot(q,[p p1 p2])
legend('Inverse Function','1st Order Approximation','2nd Order Approximation','Location','NE')
xlim([1.7 2.3])
set(gca,'xTick',[1.7 2.0 2.3])
set(gca,'yTick',[0 1 2 3])
title('Illustration of Inverse Function Theorem')
xlabel('$q$')
ylabel('$p$')


%% Implicit Function Theorem

n = 100;
x = nodeunif(n,0,2);
F = @(y) y.^5+y.^3-x.^2-1;
y = ones(n,1);
y = broyden(F,y);
f1 = 1/4;
f2 = 3/64;
y1 = 1+f1*(x-1);
y2 = y1 + 0.5*f2*(x-1).^2;

figure
plot(x,[y y1 y2])
legend('Implicit Function','1st Order Approximation','2nd Order Approximation')
xlim([0 2])
ylim([0.7 1.3])
set(gca,'xTick',[0 1 2])
set(gca,'yTick',[0.7 1.0 1.3])
title('Illustration of Implicit Function Theorem')
xlabel('$x$')
ylabel('$y$')


%% Mean Value Theorem

a = -2; b = 1;
f = @(x) x.*sin(x);
davg = (f(b)-f(a))/(b-a);
g = @(x)  x*cos(x)+sin(x)-davg;
xbar = broyden(g,0.5);
ybar = f(xbar);
x = nodeunif(100,a,b);
figure
hold on
plot(x,f(x))
plot([a b],[f(a) f(b)],'r')
plot(x,ybar + davg*(x-xbar),'r:')
yl = ylim;
plotvdash(xbar,f(xbar))
plottext(xbar+0.05,yl(1),'$\bar x$')
plottext(-1.1,1,'$f(x)$')
plottext(a,yl(1),'$a$','center','top')
plottext(b,yl(1),'$b$','center','top')
set(gca,'XTick',[])
set(gca,'YTick',[])
title('Mean Value Theorem')
xlabel('$x$')
ylabel('$y$')


%% SAVE FIGURES
printfigures(mfilename)