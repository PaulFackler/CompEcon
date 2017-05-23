% DEMAPP04 Approximating Runge's function and other comparisons
disp('DEMAPP04 Approximating Runge''s function and other comparisons')
close all

runge=inline(' 1./(1+25*x.^2)','x');

n = 13;
a = -5;
b =  5;

% Define family for evenly spaced breakpoints
breaks  = linspace(a,b,n)';
fspace1 = fundef({'spli',breaks});
c1      = funfitf(fspace1,runge);
% Define family for unevenly spaced breakpoints (note use of fundef rather than fundefn)
breaks  = a; 
for i=1:fix(n/2)-1, breaks = [breaks;breaks(end)/2]; end
breaks  = [breaks;0;flipud(-breaks)];
fspace2 = fundef({'spli',breaks});
c2      = funfitf(fspace2,runge);

x = nodeunif(501,a,b);

figure(1)
plot(x,feval(runge,x),':',x,funeval(c1,fspace1,x),'--', ...
     x,funeval(c2,fspace2,x),'-',breaks,-.17,'kx'); 
title('Runge''s Function with Spline Approximations');
xlabel('x')
ylabel('y')
legend('Runge','Even Spacing','Uneven Spacing')

figure(2)
plot(x,feval(runge,x)-funeval(c1,fspace1,x),'--', ...
     x,100*(feval(runge,x)-funeval(c2,fspace2,x)),'-',breaks,-.77,'kx'); 
title('Approximation Errors for Runge''s Function');
xlabel('x')
ylabel('error')
legend('even spacing','uneven spacing  x 100');


% Compare cubic spline with Chebychev polynomial approximations
n       = length(breaks)+2;
fspace1 = fundef({'cheb',n,a,b});
c1      = funfitf(fspace1,runge);
fspace2 = fundef({'spli',breaks});
c2      = funfitf(fspace2,runge);

figure(3)
plot(x,feval(runge,x),':',x,funeval(c1,fspace1,x),'--',...
     x,funeval(c2,fspace2,x),'-',breaks,-.17,'kx'); 
title('Runge''s Function with Chebyshev and Spline Approximations');
xlabel('x')
ylabel('y')
legend('Runge','Chebyshev','Cubic Spline')
text(-4.5,.75,{'Approximation Orders';['Chebyshev: ' num2str(n)]; ...
    ['Spline: ' num2str(length(breaks)+2)]});

figure(4)
plot(x,feval(runge,x)-funeval(c1,fspace1,x),'--', ...
     x,100*(feval(runge,x)-funeval(c2,fspace2,x)),'-',breaks,-.77,'kx'); 
title('Approximation Errors for Runge''s Function');
xlabel('x')
ylabel('y')
legend('Chebyshev','Cubic Spline')
text(-4.5,-.4,{'Approximation Orders';['Chebyshev: ' num2str(n)]; ...
    ['Spline: ' num2str(length(breaks)+2)]});
text(1.75,-.7,'Spline errors 100x')

% Compare polynomial interpolation with Chebychev and evenly spaced nodes
xnodes  = linspace(a,b,n)';
fspace1 = fundefn('cheb',n,a,b);
c1      = funfitf(fspace1,runge);
fspace2 = fundefn('cheb',n,a,b);
c2      = funfitxy(fspace2,xnodes,feval(runge,xnodes));

figure(5)
plot(x,feval(runge,x),':',x,funeval(c1,fspace1,x),'--',x,funeval(c2,fspace2,x),'-'); 
title('Runge''s Function Using Chebyshev and Uniform Nodes');
xlabel('x')
ylabel('y')
legend('Runge','Chebyshev','Uniform')
ylim([-1 2])


% Compare polynomial approximations with different numbers of nodes
testf   = inline('exp(-0.5*x.^2)');
n       = 10;
fspace1 = fundefn('cheb',n,a,b);
c1      = funfitf(fspace1,testf);
n       = 30;
fspace2 = fundefn('cheb',n,a,b);
c2      = funfitf(fspace2,testf);
x       = linspace(a,b,401)';

figure(6)
plot(x,(feval(testf,x)-funeval(c1,fspace1,x)),':', ...
     x,10e4*(feval(testf,x)-funeval(c2,fspace2,x)),'-'); 
title(['Errors in Approximating f(x)=exp(-0.5*x^2)']);
xlabel('x')
ylabel('y')
legend('n=10','n=30')
text(2.5,-.07,'n=30 10000x')

prtfigs(mfilename,'Approximation Errors for Runge''s Function',2)
