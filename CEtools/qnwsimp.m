% QNWSIMP Computes multivariate Simpson quadrature nodes and weights
% USAGE
%   [x,w] = qnwsimp(n,a,b);
% INPUTS
%   n   : 1 by d vector of number of nodes per dimension
%   a   : 1 by d vector of left endpoints
%   b   : 1 by d vector of right endpoints
% OUTPUTS
%   x   : prod(n) by d matrix of nodes
%   w   : prod(n) by 1 vector of weights
% 
% Note: each n(i) must be an odd number
%
% To compute integral If of f over interval [a,b], write a
% function f that returns an m-vector of values when passed an
% m by d matrix, and write [x,w]=qnwsimp(n,a,b); intf=w'*f(x);
%
% Integration will be exact for polynomial of order 3 and less.
%
% USES: ckron, gridmake

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,w]=qnwsimp(n,a,b)

d = length(n);
if nargin<2 | isempty(a), a=-ones(1,d); end
if nargin<3 | isempty(b), b= ones(1,d); end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwsimp1(n(i),a(i),b(i));
end
w = ckron(w(d:-1:1));
x = gridmake(x);

return


% QNWSIMP1 Computes univariate Simpson quadrature nodes and weights
% USAGE
%    [x,w] = qnwsimp1(n,a,b)
% INPUTS
%   n   : number of nodes
%   a   : left endpoint
%   b   : right endpoint
% OUTPUTS
%   x   : n by 1 matrix of nodes
%   w   : n by 1 vector of weights
% 
% Note: n must be an odd number

function [x,w]=qnwsimp1(n,a,b)

if rem(n,2) == 0
  warning('n must be an odd integer - increasing by 1');
  n=n+1;
end
dx = (b-a)/(n-1); x = (a:dx:b)';
w = reshape([2;4]*ones(1,(n+1)/2),n+1,1);
w = w(1:n); w(1) = 1; w(n) = 1;
w = (dx/3)*w;
