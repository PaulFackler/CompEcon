% QNWTRAP Computes multivariate trapezoid rule quadrature nodes and weights
% USAGE
%   [x,w] = qnwtrap(n,a,b)
% INPUTS
%   n   : 1 by d vector of number of nodes per dimension
%   a   : 1 by d vector of left endpoints
%   b   : 1 by d vector of right endpoints
% OUTPUTS
%   x   : prod(n) by d matrix of nodes
%   w   : prod(n) by 1 vector of weights
% 
% To compute integral If of f over interval [a,b], write a
% function f that returns an m-vector of values when passed an
% m by d matrix, and write [x,w]=qnwtrap(n,a,b); intf=w'*f(x);
%
% Integration will be exact for polynomial of order 1 and less.
%
% USES: ckron, gridmake

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,w] = qnwtrap(n,a,b)

d = length(n);
if nargin<2 | isempty(a), a=-ones(1,d); end
if nargin<3 | isempty(b), b= ones(1,d); end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwtrap1(n(i),a(i),b(i));
end
w = ckron(w(d:-1:1));
x = gridmake(x);

return


% QNWTRAP1 Computes univariate trapezoid rule quadrature nodes and weights
% USAGE
%    [x,w] = qnwtrap1(n,a,b);
% INPUTS
%   n   : number of nodes
%   a   : left endpoint
%   b   : right endpoint
% OUTPUTS
%   x   : n by 1 matrix of nodes
%   w   : n by 1 vector of weights

function [x,w]=qnwtrap1(n,a,b)

if n < 1
  error('n must be greater than one');
end
dx = (b-a)/(n-1);
x = (a:dx:b)';
w = dx*ones(n,1);
w([1;n]) = 0.5*w([1;n]);
