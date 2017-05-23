% DEMQUA05 Compares Chebyshev and Legendre quadrature nodes and weights
function demqua05
close all

disp(' ')
disp('DEMQUA05 Compares Chebyshev and Legendre quadrature nodes and weights')
 
n=15;
a=-1;
b=1;
[x1,w1]=qnwcheb(n,a,b);
[x2,w2]=qnwlege(n,a,b);

figure(1)
plot(x1,w1,'x', x2,w2,'.')
legend('Chebyshev','Legendre')
title('Quadrature Nodes and Weights')
xlabel('x')
ylabel('w') 

prtfigs(mfilename)