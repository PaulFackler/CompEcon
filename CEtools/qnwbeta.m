% QNWBETA Computes quadrature nodes and weights for Beta(a,b) distribution
% USAGE
%   [x,w] = qnwbeta(n,a,b);
% INPUTS
%   n   : number of quadrature points desired
%   a,b : Beta distribution parameters  
% OUTPUTS
%   x   : prod(n) by d matrix of nodes
%   w   : prod(n) by 1 vector of weights
%
% For multivariate problems with independent Beta random variables
% pass n, a, and b as 1 by d vectors
%
% To compute expectation of f(x), where x is Beta(a,b), write a
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=qnwbeta(n,a,b); E[f]=w'*f(x);
%
% Default parameter values are a=b=1 (uniform distribution)
%
% USES: ckron, gridmake

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu
 
function [x,w] = qnwbeta(n,a,b)

d = length(n);
if nargin<2 a=ones(1,d); end
if nargin<3 b=ones(1,d); end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwbeta1(n(i),a(i),b(i));
end
w = ckron(w(d:-1:1));  % use reverse ordered tensor product
x = gridmake(x);

return


% QNWBETA1  Computes quadrature nodes and weights for Beta(a,b) distribution
% USAGE
%   [x,w] = qnwbeta1(n,a,b);
% INPUTS
%   n   : number of quadrature points desired
%   a,b : Beta distribution parameters
% OUTPUTS
%   x   : quadrature nodes
%   w   : quadrature weights

% Based on an algorithm in W.H. Press, S.A. Teukolsky, W.T. Vetterling
% and B.P. Flannery, "Numerical Recipes in FORTRAN", 2nd ed.  Cambridge
% University Press, 1992.

function [x,w] = qnwbeta1(n,a,b)

a=a-1;
b=b-1;
maxit=25;
x=zeros(n,1);
w=zeros(n,1);
for i=1:n
  % Reasonable starting values 
  switch i
  case 1
    an=a/n;
    bn=b/n;
    r1=(1+a)*(2.78/(4+n*n)+0.768*an/n);
    r2=1+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
    z=1-r1/r2;
  case 2 
    r1=(4.1+a)/((1+a)*(1+0.156*a));
    r2=1+0.06*(n-8)*(1+0.12*a)/n;
    r3=1+0.012*b*(1+0.25*abs(a))/n;
    z=z-(1-z)*r1*r2*r3;
  case 3
    r1=(1.67+0.28*a)/(1+0.37*a);
    r2=1+0.22*(n-8)/n;
    r3=1+8*b/((6.28+b)*n*n);
    z=z-(x(1)-z)*r1*r2*r3;
  case n-1
    r1=(1+0.235*b)/(0.766+0.119*b);
    r2=1/(1+0.639*(n-4)/(1+0.71*(n-4)));
    r3=1/(1+20*a/((7.5+a)*n*n));
    z=z+(z-x(n-3))*r1*r2*r3;
  case n
    r1=(1+0.37*b)/(1.67+0.28*b);
    r2=1/(1+0.22*(n-8)/n);
    r3=1/(1+8*a/((6.28+a)*n*n));
    z=z+(z-x(n-2))*r1*r2*r3;
  otherwise
    z=3*x(i-1)-3*x(i-2)+x(i-3);
  end
  ab=a+b;  
  % root finding iterations 
  for its=1:maxit
    temp=2+ab;
    p1=(a-b+temp*z)/2;
    p2=1;
    for j=2:n
      p3=p2;
      p2=p1;
      temp=2*j+ab;
      aa=2*j*(j+ab)*(temp-2);
      bb=(temp-1)*(a*a-b*b+temp*(temp-2)*z);
      c=2*(j-1+a)*(j-1+b)*temp;
      p1=(bb*p2-c*p3)/aa;
    end
    pp=(n*(a-b-temp*z)*p1+2*(n+a)*(n+b)*p2)/(temp*(1-z*z));
    z1=z;
    z=z1-p1./pp;
    if abs(z-z1)<3e-14 break; end
  end
  if its>=maxit
    error('failure to converge in qnwbeta1')
  end
  x(i)=z;
  w(i)=temp/(pp*p2);
end
x=(1-x)/2;
w=w*exp(gammaln(a+n)+gammaln(b+n)-gammaln(n+1)-gammaln(n+ab+1));
w=w/(2*exp(gammaln(a+1)+gammaln(b+1)-gammaln(ab+2)));
