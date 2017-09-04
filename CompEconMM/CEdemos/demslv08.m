%% DEMSLV08 Nonlinear Complementarity Problem
%
% Solve nonlinear complementarity problem on R^2 using various methods,
% with randomly generated data.
 
function demslv08

% Preliminary tasks
demosetup(mfilename)

% Generate problem test data
z = randn(2,2)*2;
a = 1+min(z,[],2);
b = 1+max(z,[],2);
xinit = randn(2,1);

% Check Jacobian
[error,i,j] = checkjac(@f,xinit);

% Solve by applying Newton method to minmax formulation
optset('ncpsolve','type','minmax')
optset('ncpsolve','maxit',1500)
tic
[xm,zm] = ncpsolve(@f,a,b,xinit);
tm = toc;

% Solve by applying Newton method to semismooth formulation
optset('ncpsolve','type','ssmooth')
optset('ncpsolve','maxit',1500)
tic
[xs,zs] = ncpsolve(@f,a,b,xinit);
ts = toc;

% Print table header
fprintf('\n\nHundreds of seconds required to solve nonlinear complementarity\n')
fprintf('problem on R^2 using minmax and semismooth formulations, with\n')
fprintf('randomly generated bounds a = %4.2f %4.2f and b = %4.2f %4.2f\n',a,b)
fprintf('Algorithm           Time      Norm        x1     x2\n')
fprintf('Newton minmax     %6.2f  %8.0e     %5.2f  %5.2f\n',100*tm,norm(minmax(xm,a,b,zm)),xm)
fprintf('Newton semismooth %6.2f  %8.0e     %5.2f  %5.2f\n',100*ts,norm(minmax(xs,a,b,zs)),xs)


%% Function file
function [f,J] = f(x)
f = [200*x(1)*(x(2)-x(1)^2)+1-x(1); 100*(x(1)^2-x(2))];
J = [200*(x(2)-x(1)^2)-400*x(1)^2-1 200*x(1); 200*x(1) -100];