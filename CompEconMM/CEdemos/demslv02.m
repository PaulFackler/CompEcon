%% DEMSLV02 Root of Rosencrantz function
%
% Compute root of f(x1,x2)= [200*x1*(x2-x1^2)+1-x1;100*(x1^2-x2)] using
% Newton and Broyden methods with random initial value; true root is x1=1
% x2=1.

function demslv02

% Preliminary tasks
demosetup(mfilename)
  
% Randomly generate starting point
xinit = randn(2,1);

% Compute root using Newton method
optset('newton','showiters',0)
optset('newton','maxit',1500)
tic;
xn = newton(@f,xinit);
tn = toc;

% Compute root using Broyden method
optset('broyden','showiters',0)
optset('broyden','maxit',1500)
tic;
xb = broyden(@f,xinit);
tb = toc;

% Print table header
fprintf('Hundreds of seconds required to compute root of Rosencrantz function\n')
fprintf('f(x1,x2)= [200*x1*(x2-x1^2)+1-x1;100*(x1^2-x2)] via Newton and Broyden\n')
fprintf('methods, starting at x1=%4.2f x2=%4.2f\n',xinit)
fprintf('Method      Time    ||f(x)||        x1     x2\n')
fprintf('Newton  %8.2f    %8.0e     %5.2f  %5.2f\n',100*tn,norm(f(xn)),xn)
fprintf('Broyden %8.2f    %8.0e     %5.2f  %5.2f\n',100*tb,norm(f(xb)),xb)


%% Rosencrantz Function
function [fval,fjac] = f(x)
fval = [200*x(1)*(x(2)-x(1)^2)+1-x(1); 100*(x(1)^2-x(2))];
fjac = [200*(x(2)-x(1)^2)-400*x(1)^2-1 200*x(1); 200*x(1) -100];