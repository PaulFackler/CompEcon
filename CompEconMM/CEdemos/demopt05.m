%% DEMOPT05 Optimization with qnewton

% Preliminary tasks
demosetup(mfilename)

% Maximize univariate function
f = @(x) x^3-12*x^2+36*x+8;
x0 = 4;
x = qnewton(f,x0);
J = fdjac(f,x);
H = fdhess(f,x);
fprintf('Univariate Optimization\n')
fprintf('Maximum                  %8.4f\n',x)
fprintf('First Derivatve          %8.4f\n',J)
fprintf('Hecond Derivative        %8.4f\n',H)

% Maximize biivariate function
f = @(x) 5-4*x(1)^2-2*x(2)^2-4*x(1)*x(2)-2*x(2);
x0 = [-1;1];
x  = qnewton(f,x0);
J = fdjac(f,x);
E = eig(fdhess(f,x));
fprintf('\n\nBivariate Optimization\n')
fprintf('Maximum                  %8.4f   %8.4f\n',x)
fprintf('Jacobian                 %8.4f   %8.4f\n',J)
fprintf('Hessian Eigenvalues      %8.4f   %8.4f\n',E)