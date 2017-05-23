% DEMSLV04 Graphical demonstration of computing fixed point y=x^0.5 
%    using both function iteration and Newton methods
function demslv04
disp('DEMSLV04 Graphical demonstration of computing fixed point y=x^0.5')
close all

maxit = 15;
tol = 1.e-7;

alpha = 0.5;
xinit = 0.4;

x = xinit;
xf1 = []; xf2 = [];
for it=1:maxit
  xf1 = [xf1 x x];
  xf2 = [xf2 x];
  x = x^alpha;
  xf2 = [xf2 x];
end

x = xinit;
xn1 = []; xn2 = [];
for it=1:maxit
  xn1 = [xn1 x];
  xn2 = [xn2 x^alpha];
  x = x - (x-x^alpha)/(1-alpha*x^(alpha-1));
  xn1 = [xn1 x];
  xn2 = [xn2 x];
end

xplot = linspace(min([xf1 xn1])-.01,max([xf1 xn1])+.01);
yplot = xplot.^alpha;

figure(1); plot(xf1,xf2,'k',xplot,xplot,'k',xplot,yplot,'k');
title('Computing Fixed-Point via Function Iteration')

figure(2); plot(xn1,xn2,'k',xplot,xplot,'k',xplot,yplot,'k');
title('Computing Fixed-Point via Newton Iteration')

prtfigs(mfilename)
