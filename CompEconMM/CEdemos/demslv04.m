%% DEMSLV04 Fixedpoint of g(x1,x2)= [x1^2+x2^3;x1*x2-0.5]
%
% Compute fixedpoint of g(x1,x2)= [x1^2+x2^3;x1*x2-0.5] using Newton,
% Broyden, and function iteration methods with random initial value; true
% fixedpoint is x1=-0.09 x2=-0.46.  Some algorithms may fail to converge,
% depending on the initial value.

function demslv04

% Preliminary tasks
demosetup(mfilename)
warning off

% Randomly generate starting point
xinit = randn(2,1);

% Compute fixed-point using Newton method
optset('newton','showiters',0)
optset('newton','maxit',1500)
tic;
xn = newton(@f,xinit);
tn = toc;

% Compute fixed-point using Broyden method
optset('broyden','showiters',0)
optset('broyden','maxit',1500)
tic;
xb = broyden(@f,xinit);
tb = toc;

% Compute fixed-point using function iteration
optset('fixpoint','showiters',0)
optset('fixpoint','maxit',1500)
tic;
xf = fixpoint(@g,xinit);
tf = toc;

% Print table
fprintf('Hundreds of seconds required to compute fixed-point of g(x1,x2)=[x1^2+x2^3;x1*x2-0.5]\n')
fprintf('using Newton, Broyden, and function iteration methods, starting at x1=%4.2f x2=%4.2f\n',xinit)
fprintf('Method      Time  ||x-g(x)||        x1     x2\n')
fprintf('Newton  %8.2f    %8.0e     %5.2f  %5.2f\n',100*tn,norm(f(xn)),xn)
fprintf('Broyden %8.2f    %8.0e     %5.2f  %5.2f\n',100*tb,norm(f(xb)),xb)
fprintf('Function%8.2f    %8.0e     %5.2f  %5.2f\n',100*tf,norm(f(xf)),xf)


%% Function
function gval = g(x)
gval = [x(1)^2 + x(2)^3; x(1)*x(2) - 0.5];


%% Equivalent Rootfinding Formulation
function [fval,fjac] = f(x)
fval = [x(1) - x(1)^2 - x(2)^3; x(2) - x(1)*x(2) + 0.5];
fjac = [1-2*x(1) -3*x(2)^2; -x(2) 1-x(1)];