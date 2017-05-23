% QNWLEGE Computes multivariate Guass-Legendre quadrature nodes and weights
% USAGE
%    [x,w] = qnwlege(n,a,b);
% INPUTS
%   n   : 1 by d vector of number of nodes per dimension
%   a   : 1 by d vector of left endpoints
%   b   : 1 by d vector of right endpoints
% OUPUTS
%   x   : prod(n) by d matrix of nodes
%   w   : prod(n) by 1 vector of weights
% 
% To compute integral of f(x) over interval [a,b], write a
% function f that returns an m-vector of values when passed an
% m by d matrix, and write [x,w]=qnwlege(n,a,b); intf=w'*f(x);
%
% Integration will be exact for polynomial of order less than 2n.
% If a and b not specified the standard interval [-1,1] is assumed.
%
% USES: ckron, gridmake

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,w] = qnwlege(n,a,b)

d = length(n);
if nargin<2 | isempty(a), a=-ones(1,d); end
if nargin<3 | isempty(b), b= ones(1,d); end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwlege1(n(i),a(i),b(i));
end
w = ckron(w(d:-1:1));  % use reverse ordered tensor product
x = gridmake(x);

return


% QNWLEGE1 Computes univariate Gauss-Legendre quadrature nodes and weights
% USAGE
%    [x,w] = qnwlege1(n,a,b)
% INPUTS
%   n   : number of nodes
%   a   : left endpoint
%   b   : right endpoint
% OUTPUTS
%   x   : n by 1 vector of nodes
%   w   : n by 1 vector of weights

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
   if abs(z-z1)<1e-14 break; end
end
if its==maxit
   error('Maximum iterations in qnwlege1')
end
x(i) = xm-xl*z;
x(n+1-i) = xm+xl*z;
w(i) = 2*xl./((1-z.*z).*pp.*pp);
w(n+1-i) = w(i);