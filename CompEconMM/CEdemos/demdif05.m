%% DEMDIF05 Finite-Difference Jacobians and Hessians

format long

disp('Finite difference Jacobian')
f1 = @(x) [exp(x(1))-x(2);x(1)+x(2)^2;(1-x(1))*log(x(2))];
disp(fdjac(f1,[0;1]))

disp('Finite difference Hessian')
f2 = @(x) x(1)^2*exp(-x(2));
disp(fdhess(f2,[1;0]))