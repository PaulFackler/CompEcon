%% DEMDIF04 Demonstrates fdhess

% Preliminary tasks
demosetup(mfilename)

f = @(x) x(2)*exp(x(1))+x(3)^2;
x = [0;1;1];

hessanalytic  = [1 1 0;1 0 0;0 0 2];
hessnumerical = fdhess(f,x);

format long
disp('Analytic Derivative')
disp(hessanalytic)
disp('Numerical Derivative')
disp(hessnumerical)

[error,i]   = max(abs(hessnumerical-hessanalytic));
[error,j]   = max(error);
i = i(j);

fprintf('Comparing Analytic and Finite Difference Hessian\n')
fprintf('  Maximum Discrepancy %12.6e\n',error)
fprintf('  Row     %7i\n',i)
fprintf('  Column  %7i\n',j)