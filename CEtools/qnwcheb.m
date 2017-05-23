% QNWCHEB Computes multivariate Guass-Chebyshev quadrature nodes and weights
% USAGE
%    [x,w] = qnwcheb(n,a,b);
% INPUTS
%   n   : 1 by d vector of number of nodes per dimension
%   a   : 1 by d vector of left endpoints
%   b   : 1 by d vector of right endpoints
% OUTPUTS
%   x   : prod(n) by d matrix of nodes
%   w   : prod(n) by 1 vector of weights
% 
% To compute integral of f(x) over interval [a,b], write a
% function f that returns an m-vector of values when passed an
% m by d matrix, and write [x,w]=qnwcheb(n,a,b); intf=w'*f(x);
%
% USES: ckron, gridmake

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,w] = qnwcheb(n,a,b)

d = length(n);
if nargin<2 | isempty(a), a=-ones(1,d); end
if nargin<3 | isempty(b), b= ones(1,d); end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwcheb1(n(i),a(i),b(i));
end
w = ckron(w(d:-1:1));  % use reverse ordered tensor product
x = gridmake(x);


% QNWCHEB1 Computes univariate Gauss-Chebyshev quadrature nodes and weights
% USAGE
%    [x,w] = qnwcheb1(n,a,b)
% INPUTS
%   n   : number of nodes
%   a   : left endpoint
%   b   : right endpoint
% OUTPUTS
%   x   : n by 1 vector of nodes
%   w   : n by 1 vector of weights

function [x,w] = qnwcheb1(n,a,b);

x=(b+a)/2-(b-a)/2*cos(pi/n*(0.5:(n-0.5))');
w=((b-a)/n)*(cos(pi/n*((1:n)'-0.5)*(0:2:n-1)) ...
       *[1;-2./((1:2:n-2).*(3:2:n))']);

