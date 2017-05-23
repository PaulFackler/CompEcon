% DEMAPP10 Approximate function inverse via collocation
disp('DEMAPP10 Approximate function inverse via collocation')
close all

fspace = fundefn('cheb',6,1,2);      % define space of approximating functions
x      = funnode(fspace);            % compute standard nodes
Phi    = funbas(fspace,x);           % get basis matrix
c      = funfitxy(fspace,x,x);       % fit initial values using identity function
c      = broyden('fapp10',c,x,Phi);  % call rootfinding routine to compute coefficients

% plot approximation residual and errors
xplot  = nodeunif(101,1,2);
figure(1)
plot(xplot,exp(funeval(c,fspace,xplot))-xplot)
title('Residual Function: exp(f(x))-x');
xlabel('x'); ylabel('r')

figure(2)
plot(xplot,log(xplot)-funeval(c,fspace,xplot))
title('Approximation Errors for exp^{-1}(x)');
xlabel('x'); ylabel('error')

prtfigs(mfilename)
