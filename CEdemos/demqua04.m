% DEMQUA04 Compare various quadrature methods
function demqua04
close all

disp(' ')
disp('DEMQUA04 Compare various quadrature methods')

type = {'N','W','H','R'};
int = zeros(4,5);
n = [3;5;11;31];

% Define 1-d integrand and bounds
f = inline('exp(-x.*x/2)');
a = -1; b =  2;

% Perform 1-d integrations
for i=1:4;
   int(i,1) = quadrect(f,n(i)^2,a,b,'lege');
   int(i,2) = quadrect(f,n(i)^2,a,b,'trap');
   int(i,3) = quadrect(f,n(i)^2,a,b,'simp');
   int(i,4) = quadrect(f,n(i)^2,a,b,'N');
   int(i,5) = quadrect(f,n(i)^2,a,b,'W');
   int(i,6) = quadrect(f,n(i)^2,a,b,'H');
   int(i,7) = quadrect(f,n(i)^2,a,b,'R');
end

% Print out 1-d integrations
fprintf('\n Integrate exp(-x^2/2) on [%1i,%1i]\n',a,b)
fprintf('         n    Legendre   Trapezoid  Simpson    Neider.    Weyl       Haber      Random\n')
for i=1:4;
   fprintf('%10i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',[n(i)^2 int(i,:)]);
end

%textable([n int],[0 8 8 8 8 8 8 8],{'n','L','T','S','N','W','H','R'})

% Define 2-d integrand and bounds
f = inline('exp(x(:,1)+x(:,2))');
a = [ 0, 0]; b = [ 1, 2];

% Perform 2-d integrations
for i=1:4;
   int(i,1) = quadrect(f,[n(i) n(i)],a,b,'lege');
   int(i,2) = quadrect(f,[n(i) n(i)],a,b,'trap');
   int(i,3) = quadrect(f,[n(i) n(i)],a,b,'simp');
   int(i,4) = quadrect(f,n(i)^2,a,b,'N');
   int(i,5) = quadrect(f,n(i)^2,a,b,'W');
   int(i,6) = quadrect(f,n(i)^2,a,b,'H');
   int(i,7) = quadrect(f,n(i)^2,a,b,'R');
end

% Print out 2-d integrations
fprintf('\n Integration of exp(x1+x2) on unit square\n')
fprintf('         n    Legendre   Trapezoid  Simpson    Neider.    Weyl       Haber      Random\n')
for i=1:4;
   fprintf('%10i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',[n(i)^2 int(i,:)]);
end

%textable([n int],[0 8 8 8 8 8 8 8],{'n','L','T','S','N','W','H','R'})
