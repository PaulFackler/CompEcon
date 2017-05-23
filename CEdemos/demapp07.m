% DEMAPP07 Approximate y=x1/exp(x2) on [0,5]x[-1,1] using spline and Chebychev approximation
disp('DEMAPP07 Approximate y=x1/exp(x2) on [0,5]x[-1,1] using spline and Chebychev approximation')
close all
clear all

% Set degree and domain of interpolation
n = [11 11];
a = [ 0 -1];
b = [ 5  1];

% Construct Chebychev interpolant
chebproj = fundefn('cheb',n,a,b);
ccheb    = funfitf(chebproj,'fapp07');

% Construct cubic spline interpolant
spliproj = fundefn('spli',n,a,b);
cspli    = funfitf(spliproj,'fapp07');

% Construct refined uniform grid for error ploting
nplot = [41 41];
[x,xcoord] = nodeunif(nplot,a,b);

% Compute actual and fitted values on grid
[y,d] = fapp07(x);                                % actual
ycheb = funeval(ccheb,chebproj,x);                % Chebychev
ycspl = funeval(cspli,spliproj,x);                % cubic spline

% Compute fitted derivatives dy/dx1 on grid
dcheb = funeval(ccheb,chebproj,x,[1 0]);          % Chebychev
dcspl = funeval(cspli,spliproj,x,[1 0]);          % cubic spline

% Plot function to be approximated
figure(1);
yplot = reshape(y,nplot(1),nplot(2));
subplot(2,2,1); hh=surf(xcoord{1},xcoord{2},yplot);
xlabel('x1'); ylabel('x2'); zlabel('y');
title('Function to be Approximated')
set(hh,'FaceColor','interp','EdgeColor','interp')

% Plot Chebychev function approximation error
error = reshape(ycheb-y,nplot(1),nplot(2));
subplot(2,2,2); hh=surf(xcoord{1},xcoord{2},error);
xlabel('x1'); ylabel('x2'); zlabel('error');
title('Chebychev Approximation Error')
set(hh,'FaceColor','interp','EdgeColor','interp')

% Plot cubic spline function approximation error
error = reshape(ycspl-y,nplot(1),nplot(2));
subplot(2,2,3); hh=surf(xcoord{1},xcoord{2},error);
xlabel('x1'); ylabel('x2'); zlabel('error');
title('Cubic Spline Approximation Error')
set(hh,'FaceColor','interp','EdgeColor','interp')

% Plot derivative dy/dx1 to be approximated
figure(2);
yplot = reshape(d(:,1),nplot(1),nplot(2));
subplot(2,2,1); hh=surf(xcoord{1},xcoord{2},yplot);
xlabel('x1'); ylabel('x2'); zlabel('dy/dx1');
title('Derivative to be Approximated')
set(hh,'FaceColor','interp','EdgeColor','interp')

% Plot Chebychev dy/dx1 approximation error
error = reshape(dcheb-d(:,1),nplot(1),nplot(2));
subplot(2,2,2); hh=surf(xcoord{1},xcoord{2},error);
xlabel('x1'); ylabel('x2'); zlabel('error');
title('Chebychev Approximation Error')
set(hh,'FaceColor','interp','EdgeColor','interp')

% Plot cubic spline dy/dx1 approximation error
error = reshape(dcspl-d(:,1),nplot(1),nplot(2));
subplot(2,2,3); hh=surf(xcoord{1},xcoord{2},error);
xlabel('x1'); ylabel('x2'); zlabel('error');
title('Cubic Spline Approximation Error')
set(hh,'FaceColor','interp','EdgeColor','interp')
