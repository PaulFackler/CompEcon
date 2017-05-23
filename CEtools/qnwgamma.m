% QNWGAMMA Quadrature nodes and weights for Gamma(a) distribution
% USAGE
%   [x,w] = qnwgamma(n,a);
% INPUTS
%   n : the number of quadrature points desired
%   a : parameter of Gamma distribution
% OUTPUTS
%   x   : prod(n) by d matrix of nodes
%   w   : prod(n) by 1 vector of weights
%
% For multivariate problems with independent Gamma random variables
% pass n and a as 1 by d vectors
%
% To compute expectation of f(x), where x is Gamma(a), write a
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=qnwgamma(n,a); E[f]=w'*f(x);
%
% Default parameter value is a=1 (exponential distribution)
%
% USES: ckron, gridmake

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu


function [x,w] = qnwgamma(n,a)

d = length(n);
if nargin<2 a=ones(1,d); end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwgamma1(n(i),a(i));
end
w = ckron(w(d:-1:1));  % use reverse ordered tensor product
x = gridmake(x);

return


% QNWGAMMA1 Quadrature nodes and weights for Gamma(a) distribution
% USAGE
%   [x,w] = qnwgamma1(n,a);
% INPUTS
%   n : the number of quadrature points desired
%   a : parameter of Gamma distribution (default=1)
% OUTPUTS
%   x : quadrature nodes
%   w : quadrature weights
% 
% To compute expectation of f(x), where x is Gamma(a), write a
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=qnwgamm(n,a); E[f]=w'*f(x);

% Based on an algorithm in W.H. Press, S.A. Teukolsky, W.T. Vetterling
% and B.P. Flannery, "Numerical Recipes in FORTRAN", 2nd ed.  Cambridge
% University Press, 1992.

function [x,w] = qnwgamma1(n,a)

if nargin<2 a=0; else a=a-1; end

maxit=10;
factor=-exp(gammaln(a+n)-gammaln(n)-gammaln(a+1));
x=zeros(n,1);
w=zeros(n,1);
for i=1:n
  % Reasonable starting values 
  if     i==1 z=(1+a)*(3+0.92*a)/(1+2.4*n+1.8*a);
  elseif i==2 z=z+(15+6.25*a)./(1+0.9*a+2.5*n);
  else
      j=i-2;
      z=z+((1+2.55*j)./(1.9*j)+1.26*j*a./(1+3.5*j))*(z-x(j))./(1+0.3*a);
  end
  % root finding iterations 
  for its=1:maxit
    p1=1;
    p2=0;
    for j=1:n
      p3=p2;
      p2=p1;
      p1=((2*j-1+a-z)*p2-(j-1+a)*p3)./j;
    end
    pp=(n*p1-(n+a)*p2)./z;
    z1=z;
    z=z1-p1./pp;
    if abs(z-z1)<3e-14 break; end
  end
  if its>=maxit
    error('Failure to converge in qnwgamma1')
  end
  x(i)=z;
  w(i)=factor/(pp*n*p2);
end
