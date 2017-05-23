% DEMQUA03 Compares quadrature methods
disp('DEMQUA03 Compare quadrature methods')
close all

type = {'N','W','H','R'};
int = zeros(3,7,3);
n = [5;11;21;51;101;401];

% Define 1-d integrand and bounds
f1 = inline('exp(-x)');
f2 = inline('1./(1+25*x.^2)');
f3 = inline('abs(x).^0.5');

a = -1; 
b =  1;

% Perform 1-d integrations
for i=1:length(n); 
   int(i,1,1) = quadrect(f1,n(i),a,b,'trap');
   int(i,2,1) = quadrect(f1,n(i),a,b,'simp');
   int(i,3,1) = quadrect(f1,n(i),a,b,'lege');
   int(i,4,1) = quadrect(f1,n(i)^2,a,b,'N');
   int(i,5,1) = quadrect(f1,n(i)^2,a,b,'W');
   int(i,6,1) = quadrect(f1,n(i)^2,a,b,'H');
   int(i,7,1) = quadrect(f1,n(i)^2,a,b,'R');
end
for i=1:length(n); 
   int(i,1,2) = quadrect(f2,n(i),a,b,'trap');
   int(i,2,2) = quadrect(f2,n(i),a,b,'simp');
   int(i,3,2) = quadrect(f2,n(i),a,b,'lege');
   int(i,4,2) = quadrect(f2,n(i)^2,a,b,'N');
   int(i,5,2) = quadrect(f2,n(i)^2,a,b,'W');
   int(i,6,2) = quadrect(f2,n(i)^2,a,b,'H');
   int(i,7,2) = quadrect(f2,n(i)^2,a,b,'R');
end
for i=1:length(n); 
   int(i,1,3) = quadrect(f3,n(i),a,b,'trap');
   int(i,2,3) = quadrect(f3,n(i),a,b,'simp');
   int(i,3,3) = quadrect(f3,n(i),a,b,'lege');
   int(i,4,3) = quadrect(f3,n(i)^2,a,b,'N');
   int(i,5,3) = quadrect(f3,n(i)^2,a,b,'W');
   int(i,6,3) = quadrect(f3,n(i)^2,a,b,'H');
   int(i,7,3) = quadrect(f3,n(i)^2,a,b,'R');
end

% Print out 1-d integrations
fprintf('\n Integrate exp(-x^2/2) on [%1i,%1i]\n',a,b)
fprintf('         n    Legendre   Trapezoid  Simpson    Neider.    Weyl       Haber      Random\n')
for i=1:length(n); 
   fprintf('%10i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',[n(i) int(i,:,1)]);
end

fprintf('\n Integrate 1/(1+25*x.^2) on [%1i,%1i]\n',a,b)
fprintf('         n    Legendre   Trapezoid  Simpson    Neider.    Weyl       Haber      Random\n')
for i=1:length(n); 
   fprintf('%10i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',[n(i) int(i,:,2)]);
end

fprintf('\n Integrate abs(x).^0.5 on [%1i,%1i]\n',a,b)
fprintf('         n    Legendre   Trapezoid  Simpson    Neider.    Weyl       Haber      Random\n')
for i=1:length(n); 
   fprintf('%10i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',[n(i) int(i,:,3)]);
end

%textable([n int],[0 8 8 8 8 8 8 8],{'n','L','T','S','N','W','H','R'})

% Define 2-d integrand and bounds
f = inline('exp(x(:,1)+x(:,2))');
a = [ 0, 0]; b = [ 1, 2];
int=zeros(4,7);

% Perform 2-d integrations
for i=1:length(n); 
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
for i=1:length(n); 
   fprintf('%10i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',[n(i)^2 int(i,:)]);
end

%textable([n int],[0 8 8 8 8 8 8 8],{'n','L','T','S','N','W','H','R'})


% Define 2-d integrand and bounds
f = inline('exp(-x(:,1).*cos(x(:,2).^2))');
a = [ -1, -1]; b = [ 1, 1];
true = quadrect(f,[400 400],a,b,'lege');
% Perform 2-d integrations
int = zeros(1,4);
fprintf('\n Integration of exp(x1+x2) on unit square\n')
fprintf('         n    Neider.    Weyl       Haber      Random\n')
for i=1:3 
   n = 10^(2+i);
   int(1) = quadrect(f,n,a,b,'N');
   int(2) = quadrect(f,n,a,b,'W');
   int(3) = quadrect(f,n,a,b,'H');
   int(4) = quadrect(f,n,a,b,'R');
   fprintf('%10i %10.5f %10.5f %10.5f %10.5f\n',[n abs(int-true)]);
end