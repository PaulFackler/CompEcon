%% DEMSLV01 Root of f(x)=exp(-x)-1
%
% Compute root of f(x)=exp(-x)-1 using Newton and secant methods with
% random initial value; true root is x=0.

function demslv01

% Preliminary tasks
demosetup(mfilename)

% Randomly generate starting point
xinit = randn;

% Compute root using Newton method
optset('newton','showiters',0)
tic;
xn = newton(@f,xinit);
tn = toc;

% Compute root using Broyden method
optset('broyden','showiters',0)
tic;
xb = broyden(@f,xinit);
tb = toc;

% Print table
fprintf('Hundreds of seconds required to compute root of exp(-x)-1,\n')
fprintf('via Newton and Broyden methods, starting at x=%4.2f.\n',xinit)
fprintf('Method      Time      |f(x)|         x\n')
fprintf('Newton  %8.2f    %8.0e     %5.2f\n',100*tn,norm(f(xn)),xn)
fprintf('Broyden %8.2f    %8.0e     %5.2f\n',100*tb,norm(f(xb)),xb)


%% Function
function [fval,fjac] = f(x)
fval = exp(-x)-1;
fjac = -exp(-x);