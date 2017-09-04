%% QNWNORM
%
%  Discretizes multivariate normal distribution with mean mu and variance var.
%
%  Usage
%    [x,w] = qnwnorm(n,mu,var)
%  Let
%    d = dimension of random variable
%    N = prod(n) number of nodes and weights
%  Input
%    n   : 1.d number of nodes per dimension
%    mu  : 1.d mean vector (default: zeros)
%    var : d.d positive definite variance matrix (default: identity matrix)
%  Output
%    x   : N.d discrete mass points
%    w   : N.1 associated probabilities
%  Note
%    To compute Ef(X) when f is real-valued and X is Normal(mu,var) on R^d, 
%    write a Matlab function f that returns m.1 vector when passed an m.d 
%    matrix, and execute [x,w]=qnwnorm(n,mu,var); Ef=w'*f(x).

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwnorm(n,mu,var)

if max(n)==1
  x = mu;
  w = 1;
else
  d = length(n);
  if nargin<2, mu  = zeros(1,d); end
  if nargin<3, var = eye(d); end
  if size(mu,1)>1, mu=mu'; end
  x = cell(1,d);
  w = cell(1,d);
  for i=1:d
    [x{i},w{i}] = qnwnorm1(n(i));
  end
  w = ckron(w(d:-1:1));
  x = gridmake(x);
  x = x*chol(var)+mu(ones(prod(n),1),:);
end


%% QNWNORM1
%
%  Computes Gaussian quadrature nodes and weights for univariate standard
%  normal distribution.
%
%  Usage
%    [x,w] = qnwnorm1(n)
%  Input
%    n   : number of nodes per dimension
%  Output
%    x   : n.1 discrete mass points
%    w   : n.1 associated probabilities
%  Note
%    Based on an algorithm in Press, Teukolsky, Vetterling, and Flannery,
%    "Numerical Recipes in FORTRAN", 2nd ed. Cambridge U. Press, 1992.

function [x,w] = qnwnorm1(n)

if n==1
  x = 1;
  w = 1;
  return
end

maxit = 100;
pim4 = 1/pi.^0.25;
m = fix((n+1)./2);
x = zeros(n,1);
w = zeros(n,1);
for i=1:m
  % Reasonable starting values
  if i==1
    z = sqrt(2*n+1)-1.85575*((2*n+1).^(-1/6));
  elseif i==2
    z = z-1.14*(n.^0.426)./z;
  elseif i==3
    z = 1.86*z+0.86*x(1);
  elseif i==4
    z = 1.91*z+0.91*x(2);
  else
    z = 2*z+x(i-2);
  end;
  % Rootfinding iterations
  it = 0;
  while it<maxit
    it = it+1;
    p1 = pim4;
    p2 = 0;
    for j=1:n
      p3 = p2;
      p2 = p1;
      p1 = z.*sqrt(2/j).*p2-sqrt((j-1)/j).*p3;
    end
    pp = sqrt(2*n).*p2;
    z1 = z;
    z  = z1-p1./pp;
    if abs(z-z1)<1e-14; break; end;
  end;
  if it>=maxit
    error('In qnwnorm: failure to converge.')
  end
  x(n+1-i) = z;
  x(i) = -z;
  w(i) = 2./(pp.*pp);
  w(n+1-i) = w(i);
end;
w = w./sqrt(pi);
x = x*sqrt(2);