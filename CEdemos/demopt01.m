% DEMOPT01 Illustrates maximization via golden search
function demopt01
close all

disp(' ')
disp('DEMOPT01 Illustrates maximization via golden search')

x = nodeunif(500,0,3);

f = inline('x.*cos(x.^2)');

x1 = golden(f,0,1);
x2 = golden(f,2,3);
xx = golden(f,0,3);

close all
figure(1)
plot(x,f(x),'k',x1,f(x1),'k*',x2,f(x2),'k*')
title('Maximization of x cos(x^2) via golden search') 

prtfigs(mfilename,'Maximization of x cos(x^2) via golden search',1)
