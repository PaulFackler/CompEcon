%% DEMAPP07 Compute Implicit Function via Collocation

% Preliminary tasks
demosetup(mfilename)

% Residual function
resid = @(c,basis,x) funeval(c,basis,x).^-2 + funeval(c,basis,x).^-5 - 2*x;

% Approximation structure
n = 31;
a =  1;
b =  5;
basis = fundefn('cheb',n,a,b);      % define basis functions
xnode = funnode(basis);             % compute standard nodes

% Compute function inverse
c = [0.2;zeros(n-1,1)];             % fit initial values using identity function
c = broyden(resid,c,basis,xnode);   % call rootfinding routine to compute coefficients

% Plot setup
n = 1000;
x = nodeunif(n,a,b);
r = resid(c,basis,x);

% Plot implicit function
figure
plot(x,funeval(c,basis,x))
title('Implicit Function')
xlabel('$x$')
ylabel('$f(x)$')

% Plot residual
figure
hold on
plot(x,r)
plothdash([],0)
title('Functional Equation Residual')
xlabel('$x$')
ylabel('Residual')

%% SAVE FIGURES
printfigures(mfilename)