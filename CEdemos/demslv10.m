% DEMSLV10 Demonstrates bisection method
disp('DEMSLV10 Demonstrates bisection method')
close all

f = inline('x.^3-2');
a =  1;
b =  2;

tol = sqrt(eps);
sa = sign(feval(f,a));
sb = sign(feval(f,b));
if sa==sb, error ('f has same sign at endpoints'), end

i = 0;
x = (a+b)/2;
d = (b-a)/2;
while d>tol
   i = i+1;
   xx(i) = x;
   d = d/2;
   if sa == sign(feval(f,x))
      x = x+d;
   else
      x = x-d;
   end 
end

xplot=nodeunif(100,a,b);
yplot=feval(f,xplot);

close all
plot(xplot,yplot,xplot,zeros(size(xplot)),'-',xx,zeros(size(xx)),'*')
title('Computing Cube Root of 2 by Bisection','FontSize',14);

prtfigs(mfilename)