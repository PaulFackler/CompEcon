%% DEMDIF03 - Demonstrates fdjac and checkjac

function demdif03

% Preliminary tasks
demosetup(mfilename)

x = [1;1];

disp('Example 1: Function f - Analytic Jacobian Correctly Coded')
disp(' ')
[~,J] = f(x);
Jfd   = fdjac(@f,x);
disp('Analytic Derivative')
disp(J)
disp('Numerical Derivative')
disp(Jfd)
checkjac(@f,x);

disp(' ')
disp(' ')

disp('Example 2: Function ferr - Analytic Jacobian Incorrectly Coded')
disp(' ')
[~,J] = ferr(x);
Jfd   = fdjac(@ferr,x);
disp('Analytic Derivative')
disp(J)
disp('Numerical Derivative')
disp(Jfd)
checkjac(@ferr,x);

function [fval,J] = f(x)
fval = [x(2)*exp(x(1)); x(1)*cos(x(2));x(1)+x(2)^2];
J = [x(2)*exp(x(1)) exp(x(1)); cos(x(2)) -x(1)*sin(x(2)); 1 2*x(2)];

function [fval,J] = ferr(x)
fval = [x(2)*exp(x(1)); x(1)*cos(x(2));x(1)+x(2)^2];
J = [x(2)*exp(x(1)) exp(x(1)); cos(x(2)) x(1)*sin(x(2)); 1 2*x(2)];