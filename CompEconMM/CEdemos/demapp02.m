%% DEMAPP02 Approximating Functions on R^2

% This file illustrates the use of CompEcon Toolbox routines to construct
% and operate with an approximant for a function defined on a rectangle in 
% R^2.

% In particular, we construct an approximant for f(x1,x2)=cos(x1)/exp(x2)
% on [-1,1]X[-1,1].  The function posseses a closed-form, which will allow
% us to measure approximation error precisely. Of course, in most practical
% applications, the function to be approximated will not possess a known
% closed-form.

% Preliminary tasks
demosetup(mfilename)

% Function to be approximated and analytic partial derivatives
f   = @(x)  cos(x(:,1))./exp(x(:,2));
d1  = @(x) -sin(x(:,1))./exp(x(:,2));
d2  = @(x) -cos(x(:,1))./exp(x(:,2));
d11 = @(x) -cos(x(:,1))./exp(x(:,2));
d12 = @(x)  sin(x(:,1))./exp(x(:,2));
d22 = @(x)  cos(x(:,1))./exp(x(:,2));

% Set degree and domain of approximation
n = [ 6 6];     
a = [ 0 0];
b = [ 1 1];

% Choose approximation basis, compute interpolation matrix and nodes
[basis,Phi,xnodes] = fundefn('cheb',n,a,b);

% Compute basis coefficients by solving interpolation equation
c = Phi\f(xnodes);

% Having computed the basis coefficients, one may now evaluate the 
% approximant at any point x using funeval:
x = [0.5 0.5];
ffit = funeval(c,basis,x);

% ... one may also evaluate the approximant's partial derivatives at x:
dfit1 = funeval(c,basis,x,[1 0]);
dfit2 = funeval(c,basis,x,[0 1]);

% ... one may also evaluate the approximant's second own partial and cross 
% partial derivatives at x:
dfit11 = funeval(c,basis,x,[2 0]);
dfit22 = funeval(c,basis,x,[0 2]);
dfit12 = funeval(c,basis,x,[1 1]);

% Compare analytic and numerical computations
fprintf('Function Values and Derivatives of cos(x_1)/exp(x_2) at x=(0.5,0.5)\n')
fprintf('                Numerical       Analytic\n')
fprintf('Function     %12.5f   %12.5f\n',[ffit,f(x)])
fprintf('Partial 1    %12.5f   %12.5f\n',[dfit1,d1(x)])
fprintf('Partial 2    %12.5f   %12.5f\n',[dfit2,d2(x)])
fprintf('Partial 11   %12.5f   %12.5f\n',[dfit11,d11(x)])
fprintf('Partial 12   %12.5f   %12.5f\n',[dfit12,d12(x)])
fprintf('Partial 22   %12.5f   %12.5f\n',[dfit22,d22(x)])

% One may evaluate the accuracy of the Chebychev polynomial approximant by
% computing the function and derivative approximation error on a highly
% refined grid of points.

% Create refined grid for ploting:
nplot = [101 101];
[x,xcoord] = nodeunif(nplot,a,b);

% Compute approximation error on grid
yact  = f(x);
yfit  = funeval(c,basis,x);
error = reshape(yfit-yact,nplot);

% Plot approximation error
figure
mesh(xcoord{1},xcoord{2},error)
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('Error')
title('Chebychev Approximation Error')

% The plot indicates that an order 6 by 6 Chebychev approximation scheme 
% produces approximation errors no bigger in magnitude than 21x10^-6. 

% Let us now evaluate the accuracy of a cubic spline approximant.

% Set degree and domain of interpolation
n = [101 101];

% Compute Interpolation matrix and nodes
[basis,Phi,xnodes] = fundefn('spli',n,a,b);

% Compute basis coefficients
c = Phi\f(xnodes);

% Compute approximation error on grid
yact  = f(x);
yfit  = funeval(c,basis,x);
error = reshape(yfit-yact,nplot);

% Plot approximation error
figure
mesh(xcoord{1},xcoord{2},error)
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('Error')
title('Cubic Spline Approximation Error')

% The plot indicates that an 11 by 11 degree cubic spline approximation
% scheme produces approximation errors no bigger in magnitude than 1x10^-6.

%% SAVE FIGURES
printfigures(mfilename)