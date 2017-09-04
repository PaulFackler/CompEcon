%% DEMAPP01 Approximating Functions on R

% This file illustrates the use of CompEcon Toolbox routines to construct
% and operate with an approximant for a function defined on an interval of 
% the real line.

% In particular, we construct an approximant for f(x)=exp(-x) on the
% interval [-1,1].  The function possesses a closed-form, which will allow
% us to measure approximation error precisely. Of course, in most practical
% applications, the function to be approximated will not possess a known
% closed-form and will be defined only implicittly.

% Preliminary tasks
demosetup(mfilename)

% Function to be approximated and derivatives
f  = @(x)  exp(-x);
d1 = @(x) -exp(-x);
d2 = @(x)  exp(-x);

% Set degree and domain of interpolation
n = 10; 
a = -1; 
b =  1; 

% Choose approximation basis, compute interpolation matrix and nodes
[basis,Phi,xnodes] = fundefn('cheb',n,a,b);

% Compute basis coefficients by solving interpolation equation
c = Phi\f(xnodes);

% Having computed the basis coefficients, one may now evaluate the 
% approximant at x=0 using funeval:
x = 0;
ffit = funeval(c,basis,x);

% ... one may also evaluate the approximant's first and second derivatives 
% at x=0:
dfit1 = funeval(c,basis,x,1);
dfit2 = funeval(c,basis,x,2);

% ... and one may also evaluate the approximant's definite integral between 
% the lower bound a=-1 and x=0:
intfit = funeval(c,basis,x,-1);

% Compare analytic and numerical computations
fprintf('Function Values, Derivatives and Definite Integral of exp(-x) at x=0\n')
fprintf('                     Numerical       Analytic\n')
fprintf('Function          %12.8f   %12.8f\n',[ffit   f(x)])
fprintf('First Derivative  %12.8f   %12.8f\n',[dfit1  d1(x)])
fprintf('Second Derivative %12.8f   %12.8f\n',[dfit2  d2(x)])
fprintf('Definite Integral %12.8f   %12.8f\n',[intfit exp(1)-1])

% One may evaluate the accuracy of the Chebychev polynomial approximant by
% computing the function and derivative approximation error on a highly
% refined grid of points.

% Create refined grid for ploting:
nplot = 501;
x = nodeunif(nplot,a,b);

% Compute function and derivative approximation error on grid
ffit  = funeval(c,basis,x);
dfit1 = funeval(c,basis,x,1);
dfit2 = funeval(c,basis,x,2);
ferr  = ffit-f(x);
d1err = dfit1-d1(x);
d2err = dfit2-d2(x);

% Plot function approximation error
figure
hold on
plot(x,ferr)
plothdash([],0)
set(gca,'xtick',[-1 0 1])
xlabel('$x$')
ylabel('Error')
title('Chebychev Approximation Error - Function')

% Plot first derivative approximation error
figure
hold on
plot(x,d1err)
plothdash([],0)
set(gca,'xtick',[-1 0 1])
xlabel('$x$')
ylabel('Error')
title('Chebychev Approximation Error - First Derivative')

% Plot second derivative approximation error
figure
hold on
plot(x,d2err)
plothdash([],0)
set(gca,'xtick',[-1 0 1])
xlabel('$x$')
ylabel('Error')
title('Chebychev Approximation Error - Second Derivative')

% Let us now evaluate the accuracy of a cubic spline approximant.

% Set degree of interpolation
n = 15; 

% Choose approximation basis, compute interpolation matrix and nodes
[basis,Phi,xnodes] = fundefn('spli',n,a,b);

% Compute basis coefficients by solving interpolation equation
c = Phi\f(xnodes);

% Compute function and derivative approximation error on grid
ffit = funeval(c,basis,x);
dfit1 = funeval(c,basis,x,1);
dfit2 = funeval(c,basis,x,2);
ferr  = ffit-f(x);
d1err = dfit1-d1(x);
d2err = dfit2-d2(x);

% Plot function approximation error
figure
hold on
plot(x,ferr)
plothdash([],0)
set(gca,'xtick',[-1 0 1])
xlabel('$x$')
ylabel('Error')
title('Cubic Spline Approximation Error - Function')

% Plot first derivative approximation error
figure
hold on
plot(x,d1err)
plothdash([],0)
set(gca,'xtick',[-1 0 1])
xlabel('$x$')
ylabel('Error')
title('Cubic Spline Approximation Error - First Derivative')

% Plot second derivative approximation error
figure
hold on
plot(x,d2err)
plothdash([],0)
set(gca,'xtick',[-1 0 1])
xlabel('$x$')
ylabel('Error')
title('Cubic Spline Approximation Error - Second Derivative')

% The plot indicates that an degree 15 cubic spline approximation scheme
% produces approximation errors no bigger in magnitude than 5x10^-6.

%% SAVE FIGURES
printfigures(mfilename)