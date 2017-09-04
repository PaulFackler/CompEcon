%% DEMOPT08 Constrained Optimization Using nlpsolve

% Preliminary tasks
demosetup(mfilename)

f = @(x) -x(1)^2 - (x(2)-1)^2 - 3*x(1) + 2;
g = @(x) [4*x(1)+x(2);x(1)^2+x(2)*x(1)];
b = [0.5;2.0];
x = [0;1];
[x,fval,lambda,MP,exitflag] = nlpsolve(x,f,g,b);