%% QNWTRAP
%
%  Generates trapezoid rule quadrature nodes and weights
%
%  Usage
%    [x,w] = qnwtrap(n,a,b)
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
%    vector when passed an m.n matrix, and write [x,w]=qnwtrap(n,a,b); 
%    intf=w'*f(x).
%
%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwtrap(n,a,b)

d = length(n);
if nargin<2, a = zeros(1,d); end
if nargin<3, b = ones(1,d); end
if any(a>b)
  error('In qnwtrap: upper bounds must exceed lower bounds.') 
end
if any(n<2)
  error('In qnwtrap: elements of n must be integers greater than one.')
end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   dx = (b(i)-a(i))/(n(i)-1);
   ww = dx*ones(n(i),1);
   ww([1;n(i)]) = 0.5*ww([1;n(i)]);
   w{i} = ww;
   x{i} = (a(i):dx:b(i))';
end
w = ckron(w(d:-1:1));
x = gridmake(x);