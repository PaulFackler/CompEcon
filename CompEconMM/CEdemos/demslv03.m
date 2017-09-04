%% DEMSLV03 Fixedpoint of f(x) = x^0.5
%
% Compute fixedpoint of f(x) = x^0.5 using Newton, Broyden, and function
% iteration methods with random initial value; true fixedpoint is x=1.
% Some algorithms may fail to converge, depending on initial value.

function demslv03

% Preliminary tasks
demosetup(mfilename)


% Randomly generate starting point
xinit = rand+0.5;

% Compute fixed-point using Newton method
optset('newton','showiters',0)
tic;
xn = newton(@f,xinit);
tn = toc;

% Compute fixed-point using Broyden method
optset('broyden','showiters',0)
tic;
xb = broyden(@f,xinit);
tb = toc;

% Compute fixed-point using function iteration
optset('fixpoint','showiters',0)
tic;
xf = fixpoint(@g,xinit);
tf = toc;

% Print table
fprintf('Hundreds of seconds required to compute fixed-point of g(x)=sqrt(x)\n') 
fprintf('using Newton, Broyden, and function iteration methods, starting at\n')
fprintf('x=%4.2f\n',xinit)
fprintf('Method      Time      |f(x)|         x\n')
fprintf('Newton  %8.2f    %8.0e     %5.2f\n',100*tn,norm(f(xn)),xn)
fprintf('Broyden %8.2f    %8.0e     %5.2f\n',100*tb,norm(f(xb)),xb)
fprintf('Function%8.2f    %8.0e     %5.2f\n',100*tf,norm(f(xf)),xf)


%% Function
function gval = g(x)
gval = sqrt(x);


%% Equivalent Rootfinding Formulation
function [fval,fjac] = f(x)
fval = x-sqrt(x);
fjac = 1-0.5./sqrt(x);