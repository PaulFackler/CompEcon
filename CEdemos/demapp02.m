% DEMAPP02 Compare polynomial and spline approximation of function exp(-x)
disp('DEMAPP02 Compare polynomial and spline approximation of function exp(-x)')
close all

warning off

% Set domain of interpolation
a = -5;
b =  5;

% Define function
f = inline('exp(-x)');

% Construct refined uniform grid for ploting
nplot = 501;
xplot = nodeunif(nplot,a,b);

gphunif = [];
gphcheb = [];
gphcspl = [];
gphlspl = [];

% Print output header
fprintf('          Function Approximation Error            \n\n')
fprintf('Nodes     Chebychev    Uniform      Cub Spline   Lin Spline\n\n')

for n=5:5:50
   % Compute actual and fitted values on refined grid
   y = f(xplot);             % Actual values

   % Construct uniform monomial interpolant
   x = nodeunif(n,a,b);             % interpolation nodes
   c = polyfit(x,f(x),n-1);
   yfit = polyval(c,xplot);
   errunif = yfit-y;
   gphunif= [gphunif errunif];
   nrmunif = norm(errunif,inf);

   % Construct Chebychev interpolant
   S = fundefn('cheb',n,a,b);        % projection space
   c = funfitf(S,f);       % basis coefficients
   yfit = funeval(c,S,xplot);       % Chebychev interpolant
   errcheb = yfit-y;
   gphcheb= [gphcheb errcheb];
   nrmcheb = norm(errcheb,inf);

   % Construct cubic spline interpolant
   S = fundefn('spli',n,a,b);        % projection space
   c = funfitf(S,f);       % basis coefficients
   yfit = funeval(c,S,xplot);       % cubic spline interpolant
   errcspl = yfit-y;
   gphcspl= [gphcspl errcspl];
   nrmcspl = norm(errcspl,inf);

   % Construct linear spline interpolant
   S = fundefn('spli',n,a,b,1);     % projection space
   c = funfitf(S,f);                % basis coefficients
   yfit = funeval(c,S,xplot);       % linear spline interpolant
   errlspl = yfit-y;
   gphlspl= [gphlspl errlspl];
   nrmlspl = norm(errlspl,inf);

   % Print approximation error norm
   fprintf('%5i%13.1e%13.1e%13.1e%13.1e\n',n,nrmcheb,nrmunif,nrmcspl,nrmlspl)
end


% Print output header
fprintf('\n\n\n')
fprintf('                 Basis Condition Numbers         \n\n')
fprintf('Nodes     Chebychev    Uniform      Cub Spline   Lin Spline\n\n')

for n=5:5:50
   % Construct uniform monomial interpolant
   x = nodeunif(n,a,b);             % interpolation nodes
   B = ones(n,1);                   % basis matrix
   for i=1:n-1
      B = [B x.^i];
   end
   conunif = cond(B,2);

   % Construct Chebychev interpolant
   S = fundefn('cheb',n,a,b);        % projection space
   B = funbas(S);                  % basis matrix
   concheb = cond(B,2);

   % Construct cubic spline interpolant
   S = fundefn('spli',n,a,b);        % projection space
   B = funbas(S);                  % basis matrix
   concspl = cond(B,2);

   % Construct linear spline interpolant
   S = fundefn('spli',n,a,b,1);      % projection space
   B = funbas(S);                  % basis matrix
   conlspl = cond(B,2);

   % Print basis condition number
   fprintf('%5i%13.1e%13.1e%13.1e%13.1e\n',n,concheb,conunif,concspl,conlspl)
end

% Plot function and Chebychev and uniform node errors
close all
figure(1);
plot(xplot,gphcheb(:,2),xplot,gphunif(:,2));
xlabel('x'); ylabel('y');
title('Approximation Error for exp(-x)');
legend('Chebychev Nodes','Uniform Nodes')
set(gca,'Xtick',[-5 0 5]);

prtfigs(mfilename,'Approximation Error for exp(-x)',1);
