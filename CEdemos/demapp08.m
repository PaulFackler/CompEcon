% DEMAPP08 Approximate y=(exp(x1*x2)+x3^3)/100 on [0,1]^3 via spline and Chebychev approximation
disp('DEMAPP08 Approximate y=(exp(x1*x2)+x3^3)/100 on [0,1]^3 via spline and Chebychev approximation')
close all

% Set degree and domain of interpolation
n = [5 5 5];
a = [0 0 0];
b = [1 1 1];

% Construct Chebychev interpolant
chebproj = fundefn('cheb',n,a,b);
ccheb    = funfitf(chebproj,'fapp08');

% Construct cubic spline interpolant
spliproj = fundefn('spli',n,a,b);
cspli    = funfitf(spliproj,'fapp08');

% Construct refined uniform grid for error checking
ncheck = [20 20 20];
[x,xcoord] = nodeunif(ncheck,a,b);

% Compute actual and fitted values on grid
[y,d] = fapp08(x);                                % actual
ycheb = funeval(ccheb,chebproj,x);                % Chebychev
ycspl = funeval(cspli,spliproj,x);                % cubic spline

% Compute  fitted derivatives on grid
dcheb = funeval(ccheb,chebproj,x,eye(3));         % Chebychev
dcspl = funeval(cspli,spliproj,x,eye(3));         % cubic spline

% Estimate approximation error
fprintf('\nError         Chebychev   Cubic Spl')
errcheb = norm(ycheb-y,inf);
errcspl = norm(ycspl-y,inf);
fprintf('\nFunction     %10.2e  %10.2e',errcheb,errcspl)
errcheb = norm(dcheb(:,1)-d(:,1),inf);
errcspl = norm(dcspl(:,1)-d(:,1),inf);
fprintf('\ndf/x1        %10.2e  %10.2e',errcheb,errcspl)
errcheb = norm(dcheb(:,2)-d(:,2),inf);
errcspl = norm(dcspl(:,2)-d(:,2),inf);
fprintf('\ndf/x2        %10.2e  %10.2e',errcheb,errcspl)
errcheb = norm(dcheb(:,3)-d(:,3),inf);
errcspl = norm(dcspl(:,3)-d(:,3),inf);
fprintf('\ndf/x3        %10.2e  %10.2e\n',errcheb,errcspl)
