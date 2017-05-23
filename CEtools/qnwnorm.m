% QNWNORM Computes nodes and weights for multivariate normal distribution
% USAGE
%   [x,w] = qnwnorm(n,mu,var);
% INPUTS
%   n   : 1 by d vector of number of nodes for each variable
%   mu  : 1 by d mean vector
%   var : d by d positive definite covariance matrix
% OUTPUTS
%   x   : prod(n) by d matrix of evaluation nodes
%   w   : prod(n) by 1 vector of probabilities
% 
% To compute expectation of f(x), where x is N(mu,var), write a
% function f that returns m-vector of values when passed an m by d
% matrix, and write [x,w]=qnwnorm(n,mu,var); E[f]=w'*f(x);
%
% Options (Use OPTSET)
%   usesqrtm : 0/1 if 1 uses sqrtm to factorize var rather than chol
%                sqrtm produces a symmetric set of nodes that are
%                invariant to reordering.
%
% USES: ckron, gridmake

% Changes: 
%   4/21/2010 added usesqrtm option

% Copyright (c) 1997-2010, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,w] = qnwnorm(n,mu,var)
usesqrtm = optget('qnwnorm','usesqrtm',0);

d = length(n);
if nargin<3, var=eye(d); end
if nargin<2, mu=zeros(1,d); end
if size(mu,1)>1, mu=mu'; end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwnorm1(n(i));
end
w = ckron(w(d:-1:1));
x = gridmake(x);
if usesqrtm
  x = x*sqrtm(var)+mu(ones(prod(n),1),:);
else
  x = x*chol(var)+mu(ones(prod(n),1),:);
end
return


% QNWNORM1 Computes nodes and weights for the univariate standard normal distribution
% USAGE
%    [x,w] = qnwnorm1(n);
% INPUTS
%   n   : number of nodes
% OUTPUTS
%   x   : n by 1 vector of evaluation nodes
%   w   : n by 1 vector of probabilities
 
% Based on an algorithm in W.H. Press, S.A. Teukolsky, W.T. Vetterling
% and B.P. Flannery, "Numerical Recipes in FORTRAN", 2nd ed.  Cambridge
% University Press, 1992.

function [x,w] = qnwnorm1(n);

maxit = 100;
pim4 = 1/pi.^0.25;
m = fix((n+1)./2);
x = zeros(n,1);
w = zeros(n,1);
for i=1:m
   % Reasonable starting values 
   if i==1;        z = sqrt(2*n+1)-1.85575*((2*n+1).^(-1/6));
      elseif i==2; z = z-1.14*(n.^0.426)./z;
      elseif i==3; z = 1.86*z+0.86*x(1);
      elseif i==4; z = 1.91*z+0.91*x(2);
      else;        z = 2*z+x(i-2);
   end;
   % root finding iterations 
   its=0;
   while its<maxit;
      its = its+1;
      p1 = pim4;
      p2 = 0;
      for j=1:n
         p3 = p2;
         p2 = p1;
         p1 = z.*sqrt(2/j).*p2-sqrt((j-1)/j).*p3;
      end;
      pp = sqrt(2*n).*p2;
      z1 = z;
      z  = z1-p1./pp;
      if abs(z-z1)<1e-14; break; end;
   end;
   if its>=maxit
      error('failure to converge in qnwnorm1')
   end
   x(n+1-i) = z;
   x(i) = -z;
   w(i) = 2./(pp.*pp);
   w(n+1-i) = w(i);
end;
w = w./sqrt(pi);
x = x*sqrt(2);
