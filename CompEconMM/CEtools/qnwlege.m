%% QNWLEGE
%
%  Generates Guass-Legendre quadrature nodes and weights
%
%  Usage
%    [x,w] = qnwlege(n,a,b)
%  Let
%    d = dimension of integration interval
%    N = prod(n) number of nodes and weights
%  Input
%    n   : 1.d number of nodes per dimension
%    a   : 1.d lower bounds of integration interval
%    b   : 1.d upper bounds of integration interval
%  Output
%    x   : N.d quadrature nodes
%    w   : N.1 quadrature weights
%  Note
%    To compute definte integral of a real-valued function f defined on a
%    hypercube [a,b] in R^d, write a Matlab function f that returns an m.1 
%    vector when passed an m.n matrix, and write [x,w]=qnwlege(n,a,b); 
%    intf=w'*f(x).

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwlege(n,a,b)

d = length(n);
if nargin<2, a = zeros(1,d); end
if nargin<3, b = ones(1,d); end
if any(a>b)
  error('In qnwlege: upper bounds must exceed lower bounds.')
end
if any(n<2)
  error('In qnwlege: elements of n must be integers greater than one.')
end

x = cell(1,d);
w = cell(1,d);
for i=1:d
  [x{i},w{i}] = qnwlege1(n(i),a(i),b(i));
end
w = ckron(w(d:-1:1));  % use reverse ordered tensor product
x = gridmake(x);


%% QNWLEGE1 
%
%  Generates one-dimensional Guass-Legendre quadrature nodes and weights
%
%  Usage
%    [x,w] = qnwlege(n,a,b)
%  Input
%    n   : number of nodes
%    a   : lower bound of integration interval
%    b   : upper bound of integration interval
%  Output
%    x   : n.1 quadrature nodes
%    w   : n.1 quadrature weights
%  Note
%    Based on an algorithm in Press, Teukolsky, Vetterling, and Flannery,
%    "Numerical Recipes in FORTRAN", 2nd ed. Cambridge U. Press, 1992.

function [x,w] = qnwlege1(n,a,b)

maxit = 100;
m = fix((n+1)/2);
xm = 0.5*(b+a);
xl = 0.5*(b-a);
x = zeros(n,1);
w = x;
i = (1:m)';
z = cos(pi*(i-0.25)./(n+0.5));
for its=1:maxit
  p1 = 1;
  p2 = 0;
  for j=1:n
    p3 = p2;
    p2 = p1;
    p1 = ((2*j-1)*z.*p2-(j-1)*p3)./j;
  end
  pp = n*(z.*p1-p2)./(z.*z-1);
  z1 = z;
  z = z1-p1./pp;
  if abs(z-z1)<1e-14
    break;
  end
end
if its==maxit
  error('In qnwlege: failure to converge.')
end
x(i) = xm-xl*z;
x(n+1-i) = xm+xl*z;
w(i) = 2*xl./((1-z.*z).*pp.*pp);
w(n+1-i) = w(i);