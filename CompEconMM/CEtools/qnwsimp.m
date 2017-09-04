%% QNWSIMP 
%
%  Generates Simpson's rule quadrature nodes and weights
%
%  Usage
%    [x,w] = qnwsimp(n,a,b)
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
%    vector when passed an m.n matrix, and write [x,w]=qnwsimp(n,a,b); 
%    intf=w'*f(x).

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwsimp(n,a,b)

d = length(n);
if nargin<2, a = zeros(1,d); end
if nargin<3, b = ones(1,d); end
if any(a>b)
  error('In qnwsimp: upper bounds must exceed lower bounds.') 
end
if any(n<2)
  error('In qnwsimp: elements of n must be integers greater than one.')
end

x = cell(1,d);
w = cell(1,d);
for i=1:d
  if rem(n(i),2)==0
    warning('In qnwsimp: n must be odd integer - increasing by 1.')
    n(i) = n(i)+1;
  end
  dx = (b(i)-a(i))/(n(i)-1);
  ww = reshape([2;4]*ones(1,(n(i)+1)/2),n(i)+1,1);
  ww = ww(1:n(i));
  ww([1;n(i)]) = 1;
  w{i} =  (dx/3)*ww;
  x{i} = (a(i):dx:b(i))';
end
w = ckron(w(d:-1:1));
x = gridmake(x);


%% QNWSIMP1
%
%  Generates Simpson's rule quadrature nodes and weights for computing the
%  definite integral of a real-valued function defined on an interval [a,b]
%  in R.
%
%  Usage
%    [x,w] = qnwsimp1(n,a,b)
%  Input
%    n   : number of nodes (must be an odd positive integer)
%    a   : lower bound of integration interval
%    b   : upper bound of integration interval
%  Output
%    x   : n.1 quadrature nodes
%    w   : n.1 quadrature weights

function [x,w] = qnwsimp1(n,a,b)

dx = (b-a)/(n-1); 
x = (a:dx:b)';
w = reshape([2;4]*ones(1,(n+1)/2),n+1,1);
w = w(1:n); 
w(1) = 1; 
w(n) = 1;
w = (dx/3)*w;